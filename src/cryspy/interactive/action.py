from logging import getLogger
import os
from typing import Callable

from ..start import cryspy_init
from ..util.utility import (
    backup_cryspy,
    clean_cryspy,
    get_version,
    set_logger,
)
from .restart_interact import restart_interact


# ---------- logger
set_logger(
    logfile='log_cryspy',
    errfile='err_cryspy',
)
logger = getLogger('cryspy')


def _log_banner():
    logger.info(f'\n\n\nCrySPY {get_version()}\n\n')


def initialize():
    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise IOError('lock_cryspy file exists')
    with open('lock_cryspy', 'w'):
        pass    # create vacant file

    try:
        # ---------- check mode
        if os.path.isfile('data/db_data/rslt_data.db'):
            logger.error(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )
            raise IOError(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )

        # ---------- banner
        _log_banner()

        # ---------- initialize
        if not os.path.isfile('cryspy.stat'):
            cryspy_init.initialize()
        else:
            logger.error(
                'cryspy.stat file exists. '
                'Clean files to start from the beginning.'
            )

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


def backup():
    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise IOError('lock_cryspy file exists')
    with open('lock_cryspy', 'w'):
        pass    # create vacant file

    try:
        # ---------- check mode
        if os.path.isfile('data/db_data/rslt_data.db'):
            logger.error(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )
            raise IOError(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )

        # ---------- banner
        _log_banner()

        # ---------- backup
        backup_cryspy()

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


def clean(skip_yes=False):
    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise IOError('lock_cryspy file exists')
    with open('lock_cryspy', 'w'):
        pass    # create vacant file

    try:
        # ---------- check mode
        if os.path.isfile('data/db_data/rslt_data.db'):
            logger.error(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )
            raise IOError(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )

        # ---------- banner
        _log_banner()

        # ---------- clean
        clean_cryspy(skip_yes)

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


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

    try:
        # ---------- check mode
        if os.path.isfile('data/db_data/rslt_data.db'):
            logger.error(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )
            raise SystemExit(1)

        # ---------- banner
        _log_banner()

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

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')
