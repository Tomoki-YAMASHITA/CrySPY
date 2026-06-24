from logging import getLogger
import os
from typing import Callable

from ..start import cryspy_init
from ..util.utility import set_logger, backup_cryspy, clean_cryspy
from .restart_interact import restart_interact


# ---------- logger
set_logger(
    logfile='log_cryspy',
    errfile='err_cryspy',
)
logger = getLogger('cryspy')


def initialize():
    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise IOError('lock_cryspy file exists')
    with open('lock_cryspy', 'w'):
        pass    # create vacant file
    # ---------- initialize
    if not os.path.isfile('cryspy.stat'):
        cryspy_init.initialize()
    else:
        logger.error('cryspy.stat file exists. Clean files to start from the beginning.')
    os.remove('lock_cryspy')


def backup():
    backup_cryspy()


def clean(skip_yes=False):
    clean_cryspy(skip_yes)


def restart(
        njob: int,
        calculator: Callable,
        optimizer: str,
        symmetry: bool = True,
        fmax: float = 0.01,
        steps: int = 2000,
    ) -> None:
    """
    Restart the CrySPY process.

    Args:
        njob (int): Number of jobs to run.
        calculator (Callable): Calculator function to use.
        optimizer (str): Optimizer to use ('BFGS', 'LBFGS', 'FIRE').
        symmetry (bool, optional): Whether to use symmetry. Default is True.
        fmax (float, optional): Maximum force. Default is 0.01.
        steps (int, optional): Number of steps. Default is 2000.

    Raises:
        SystemExit: If the 'lock_cryspy' file exists.
    """

    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)

    # ---------- check args
    if njob < 0:
        logger.error('njob must be positive.')
        raise ValueError('njob must be positive.')
    if optimizer not in ['BFGS', 'LBFGS', 'FIRE']:
        logger.error(f'optimizer = {optimizer} is not supported.')
        raise ValueError(f'optimizer = {optimizer} is not supported.')

    # ---------- create lock file
    with open('lock_cryspy', 'w'):
        pass    # create vacant file

    # ---------- restart
    if os.path.isfile('cryspy.stat'):
        restart_interact(
            njob,
            calculator,
            optimizer,
            symmetry,
            fmax,
            steps,
        )
    else:
        logger.error('cryspy.stat file does not exist.')
    os.remove('lock_cryspy')
