'''
Initialize CrySPY high-throughput mode
'''

from contextlib import closing
from importlib.metadata import version
from logging import getLogger
import os

from ...IO import pkl_data, write_input
from ...IO.load_struc import load_structures, move_loaded_files
from ...IO.read_input import ReadInput
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
    init_struc_data, loaded_files = load_initial_structures_ht(rin)
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

    # ---------- composition constraints
    output_composition_constraints_ht(rin, generation=1)

    # ---------- save loaded initial structures
    save_initial_structures_ht(rin, init_struc_data)

    # ---------- move loaded files
    if loaded_files:
        move_loaded_files(loaded_files)

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


def output_composition_constraints_ht(rin, generation: int) -> None:
    """Log composition constraints and save the composition window."""
    if rin.algo != 'EA-vc':
        return

    logger.info('# ---------- Composition constraints')
    if rin.min_comp is None and rin.max_comp is None:
        logger.info('Composition constraints are not set.')
        return

    from ...util.struc_util import get_feasible_composition
    feasible_comp = get_feasible_composition(rin.min_comp, rin.max_comp)
    logger.info('Feasible composition range:')
    for symbol, (lower, upper) in zip(rin.atype, feasible_comp):
        logger.info(f'    {symbol}: {lower:.3f} - {upper:.3f}')

    if len(rin.atype) not in (2, 3):
        return

    from ...util.visual_util import save_composition_window
    save_composition_window(
        atype=rin.atype,
        gen=generation,
        min_comp=rin.min_comp,
        max_comp=rin.max_comp,
        fig_format=rin.fig_format,
        axis_order=rin.axis_order,
    )
    logger.info(
        'Composition window saved as '
        f'    ./data/convex_hull/composition_window_{generation}.{rin.fig_format}'
    )


def load_initial_structures_ht(rin):
    """Load and validate initial structures for high-throughput mode."""

    # ---------- load initial structures
    try:
        init_struc_data, loaded_files = load_structures(rin)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- return
    return init_struc_data, loaded_files


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
