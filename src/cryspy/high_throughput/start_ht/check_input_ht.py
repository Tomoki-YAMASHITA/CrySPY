from logging import getLogger
from pathlib import Path


logger = getLogger('cryspy')


def check_input_ht(rin):
    """Check input for high-throughput mode."""

    # ---------- algorithm
    if rin.algo not in ('RS', 'EA', 'EA-vc'):
        logger.error(
            'High-throughput mode currently supports only RS, EA, and EA-vc'
        )
        raise SystemExit(1)

    # ---------- structure generation
    if rin.use_find_wy:
        logger.error(
            'High-throughput mode currently does not support find_wy'
        )
        raise SystemExit(1)

    # ---------- calculation code
    if rin.calc_code != 'ASE':
        logger.error(
            'High-throughput mode currently supports only ASE'
        )
        raise SystemExit(1)

    # ---------- check point 1
    if rin.stop_chkpt == 1:
        logger.info('Stop at check point 1')
        raise SystemExit()


def check_calc_file_ht(rin):
    """Check calculation file for high-throughput mode."""

    # ---------- calculation file
    calc_path = Path('calc_in') / rin.ase_python
    if not calc_path.is_file():
        logger.error(f'{calc_path} not found')
        raise SystemExit(1)
