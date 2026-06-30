from logging import getLogger
from pathlib import Path


logger = getLogger('cryspy')


def check_input(rin):
    """Check input for high-throughput mode."""

    # ---------- algorithm
    if rin.algo != 'RS':
        logger.error(
            'High-throughput mode currently supports only RS'
        )
        raise SystemExit(1)

    # ---------- calculation code
    if rin.calc_code != 'ASE':
        logger.error(
            'High-throughput mode currently supports only ASE'
        )
        raise SystemExit(1)

    # ---------- calculation file
    calc_path = Path('calc_in') / rin.ase_python
    if not calc_path.is_file():
        logger.error(f'{calc_path} not found')
        raise SystemExit(1)

    # ---------- check point 1
    if rin.stop_chkpt == 1:
        logger.info('Stop at check point 1')
        raise SystemExit()