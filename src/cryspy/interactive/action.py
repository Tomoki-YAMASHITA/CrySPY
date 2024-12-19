from logging import getLogger
import os

from ..start import cryspy_init
from ..util.utility import set_logger, backup_cryspy, clean_cryspy

from .restart_interact import restart_interact
from .rslt_energy import display_rslt
from .view_atom import get_ase_atoms
from .energy_plot import energy_plot


set_logger(
    logfile='log_cryspy',
    errfile='err_cryspy',
)
#set_logger()
logger = getLogger('cryspy')


def initialize():
    # ---------- ignore MPI
    comm = None
    mpi_rank = 0
    mpi_size = 1
    # ---------- lock file
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise IOError('lock_cryspy file exists')
    with open('lock_cryspy', 'w'):
        pass    # create vacant file
    # ---------- initialize
    if not os.path.isfile('cryspy.stat'):
        cryspy_init.initialize(comm, mpi_rank, mpi_size)
    else:
        logger.error('cryspy.stat file exists. Clean files to start from the beginning.')
        raise IOError('cryspy.stat file exists. Clean files to start from the beginning.')
    os.remove('lock_cryspy')


def backup():
    backup_cryspy()


def clean(skip_yes=False):
    clean_cryspy(skip_yes)


def restart(
        njob: int,
        calculator,
        optimizer,
        symmetry=True,
        mask=None,
        fmax=0.01,
        steps=2000,
    ):
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    with open('lock_cryspy', 'w'):
        pass    # create vacant file
    if os.path.isfile('cryspy.stat'):
        restart_interact(
            njob,
            calculator,
            optimizer,
            symmetry,
            mask,
            fmax,
            steps,
        )
    else:
        logger.error('cryspy.stat file does not exist.')
    os.remove('lock_cryspy')


def show_rslt(cid='all'):
    display_rslt(cid)


def get_atoms(status, cid):
    return get_ase_atoms(status, cid)


def plot_e(y_max=2.0, y_min=-0.5):
    energy_plot(y_max, y_min)