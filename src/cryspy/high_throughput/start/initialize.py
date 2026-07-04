'''
Initialize CrySPY high-throughput mode
'''

from importlib.metadata import version
from logging import getLogger
import os

from ...IO import pkl_data, write_input
from ...IO.read_input import ReadInput
from ..db.sqlite import initialize_db
from ..worker.controller_rs import generate_structures_rs
from .check_input import check_input


logger = getLogger('cryspy')


def initialize():
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
    check_input(rin)

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
    initialize_db()
    logger.info('Database: data/db_data/rslt_data.db')

    # ---------- generate and insert initial structures
    logger.info('# ---------- Initial structure generation')
    generate_structures_rs(
        rin,
        range(1, rin.tot_struc + 1),
    )

    # ---------- return
    return rin
