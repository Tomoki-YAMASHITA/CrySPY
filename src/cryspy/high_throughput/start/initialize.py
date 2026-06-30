'''
Initialize CrySPY high-throughput mode
'''

from importlib.metadata import version
from logging import getLogger
import os

import numpy as np

from ...IO import write_input
from ...IO.read_input import ReadInput
from ...RS.rs_gen import gen_random
from ...util.utility import get_version
from ..db.record import insert_init_struc
from ..db.sqlite import connect_db, initialize_db
from .check_input import check_input


logger = getLogger('cryspy')


def initialize():
    # ---------- start
    logger.info('\n\n\nStart CrySPY high-throughput mode ' + get_version() + '\n\n')

    # ---------- check versions
    logger.info('# ---------- Library version info')
    logger.info(f'ase version: {version("ase")}')
    logger.info(f'pymatgen version: {version("pymatgen")}')
    logger.info(f'pyxtal version: {version("pyxtal")}')

    # ---------- read input
    logger.info('# ---------- Read input file, cryspy.in')
    try:
        rin = ReadInput()
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- check input
    check_input(rin)

    # ---------- write input
    write_input.out_input(rin)

    # ---------- make data directory
    os.makedirs('data/db_data', exist_ok=True)

    # ---------- RNG
    rng = None
    if rin.seed is not None:
        logger.info('# ---------- Initialize RNG with seed')
        rng = np.random.default_rng(rin.seed)
        logger.info(f'RNG seed: {rin.seed}')

    # ---------- initial structure generation
    logger.info('# ---------- Initial structure generation')
    init_struc_data = gen_random(
        rin=rin,
        nstruc=rin.tot_struc,
        id_offset=1,
        comm=None,
        mpi_rank=0,
        mpi_size=1,
        rng=rng,
    )

    # ---------- initialize database
    logger.info('# ---------- Initialize database')
    initialize_db()
    logger.info('Database: data/db_data/rslt_data.db')

    # ---------- insert initial structures
    logger.info('# ---------- Insert initial structures')
    with connect_db() as conn:
        for struc in init_struc_data.values():
            insert_init_struc(conn, struc, rin.nat)

    # ---------- return
    return rin