'''
Initialize CrySPY high-throughput mode
'''

from contextlib import closing
from importlib.metadata import version
from logging import getLogger
import os

from ...IO import pkl_data, write_input
from ...IO.read_input import ReadInput
from ...util.struc_validation import validate_loaded_structures
from ..db.record_db import (
    Status,
    initialize_record,
    update_init_struc,
)
from ..db.sqlite_db import connect_db, initialize_db
from ..worker.controller_ea import initialize_ea
from ..worker.controller_rs import generate_structures_rs
from .check_input_ht import check_input_ht


logger = getLogger('cryspy')


def initialize_ht():
    # ---------- check versions
    logger.info('# ---------- Library version info')
    logger.info(f'ase version: {version("ase")}')
    logger.info(f'pymatgen version: {version("pymatgen")}')
    logger.info(f'pyxtal version: {version("pyxtal")}')

    # ---------- read input
    logger.info('# ---------- Read input file, cryspy.in')
    try:
        rin = ReadInput(ht=True)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- check input
    check_input_ht(rin)

    # ---------- load initial structures
    init_struc_data = load_initial_structures_ht(rin)
    n_loaded = len(init_struc_data)

    # ---------- make data directories
    os.makedirs('data/db_data', exist_ok=True)
    os.makedirs('data/pkl_data', exist_ok=True)

    # ---------- write and save input
    write_input.out_input(rin)
    pkl_data.save_input(rin)

    # ---------- RNG
    if rin.seed is not None:
        logger.info('# ---------- Use ID-based RNG')
        logger.info(f'RNG seed: {rin.seed}')

    # ---------- initialize database
    logger.info('# ---------- Initialize database')
    initialize_db(rin.algo)
    logger.info('Database: data/db_data/rslt_data.db')

    # ---------- save loaded initial structures
    save_initial_structures_ht(rin, init_struc_data)

    # ---------- generate and insert initial structures
    logger.info('# ---------- Initial structure generation')
    if rin.algo == 'RS':
        cids = list(range(n_loaded + 1, rin.tot_struc + 1))
        generate_structures_rs(
            rin,
            cids,
        )
    elif rin.algo in ('EA', 'EA-vc'):
        cids = list(range(1, rin.n_pop + 1))
        initialize_ea(
            rin,
            cids,
        )

    # ---------- return
    return rin


def load_initial_structures_ht(rin):
    """Load and validate initial structures for high-throughput mode."""

    # ---------- load flag
    if not rin.load_struc_flag:
        # An empty pool leaves all initial structures for random generation.
        return {}

    # ---------- load initial structures
    logger.info('# ---------- Load initial structure data')
    logger.info('Load ./data/pkl_data/init_struc_data.pkl')
    try:
        init_struc_data = pkl_data.load_init_struc()
        validate_loaded_structures(
            rin,
            init_struc_data,
            id_start=1,
        )
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- return
    return init_struc_data


def save_initial_structures_ht(rin, init_struc_data):
    """Save loaded initial structures as waiting records."""

    # ---------- no loaded structures
    if not init_struc_data:
        return

    # ---------- generation
    generation = 1 if rin.algo in ('EA', 'EA-vc') else None

    # ---------- save initial structures
    with closing(connect_db()) as conn, conn:
        for cid in sorted(init_struc_data):
            initialize_record(
                conn,
                cid,
                generation=generation,
            )
            update_init_struc(
                conn,
                cid,
                init_struc_data[cid],
                rin.atype,
                rin.symprec,
                Status.WAITING,
                vc=rin.algo == 'EA-vc',
            )
