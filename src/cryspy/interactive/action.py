from logging import getLogger
import os
from typing import Callable

from ..start import cryspy_init
from ..util.utility import set_logger, backup_cryspy, clean_cryspy

from .restart_interact import restart_interact
from .manage_structure import get_ase_atoms
from .energy_plot import plot_energy, interact_plot_conv_hull, plot_conv_hull_binary, plot_conv_hull_ternary
from .show_dist import show_dist


# ---------- logger
set_logger(
    logfile='log_cryspy',
    errfile='err_cryspy',
)
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

    # ---------- crspy.stat
    with open('lock_cryspy', 'w'):
        pass    # create vacant file

    # ---------- check args
    if njob < 0:
        logger.error('njob must be positive.')
        raise ValueError('njob must be positive.')
    if optimizer not in ['BFGS', 'LBFGS', 'FIRE']:
        logger.error(f'optimizer = {optimizer} is not supported.')
        raise ValueError(f'optimizer = {optimizer} is not supported.')

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


def get_atoms(status: str, cid: int):
    return get_ase_atoms(status, cid)


def plot_E(
        title=None,
        ymax=2.0,
        ymin=-0.2,
        markersize=12,
        marker_edge_width=1.0,
        marker_edge_color='black',
        alpha=1.0,
    ):
    fig, ax = plot_energy(title, ymax, ymin, markersize, marker_edge_width, marker_edge_color, alpha)
    return fig, ax


def interactive_plot_convex_hull(cgen=None, show_unstable=0.2, ternary_style='2d'):
    interact_plot_conv_hull(cgen, show_unstable, ternary_style)


def plot_convex_hull_binary(cgen=None, show_max=0.2, label_stable=True, vmax=0.2, bottom_margin=0.02):
    fig, ax = plot_conv_hull_binary(cgen, show_max, label_stable, vmax, bottom_margin)
    return fig, ax


def plot_convex_hull_ternary(cgen=None, show_max=0.2, label_stable=True, vmax=0.2):
    fig, ax = plot_conv_hull_ternary(cgen, show_max, label_stable, vmax)
    return fig, ax


#def structure_distance(r_cut=6.0, n_max=8, l_max=2, add_str=None, add_e=None):
#    show_dist(r_cut, n_max, l_max, add_str, add_e)