from logging import getLogger
import os

from ..start import cryspy_init
from ..util.utility import set_logger, backup_cryspy, clean_cryspy

from .restart_interact import restart_interact
from .rslt_energy import display_rslt
from .view_atom import get_ase_atoms
from .energy_plot import energy_plot, energy_plot_EA, convex_hull_plot
from .show_dist import show_dist


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


def plot_e_EA(y_max=2.0, y_min=-0.5):
    energy_plot_EA(y_max, y_min)


def interactive_plot_convex_hull(show_unstable=0.05, ternary_style='2d'):
    convex_hull_plot(show_unstable, ternary_style)


def structure_distance(r_cut=6.0, n_max=8, l_max=2, add_str=None, add_e=None):
    show_dist(r_cut, n_max, l_max, add_str, add_e)