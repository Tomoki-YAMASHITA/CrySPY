#!/usr/bin/env python3
'''
Energy plot script
'''

from logging import getLogger
import os

from cryspy.IO import pkl_data
from cryspy.start import cryspy_restart
from cryspy.util.utility import set_logger
from cryspy.util.visual_util import (
    plot_energy_RS,
    plot_energy_EA,
    draw_convex_hull_binary,
    draw_convex_hull_ternary,
)


def main():
    # ---------- logger
    set_logger(
        noprint=False,
        debug=False,
        logfile='log_cryspy',
        errfile='err_cryspy',
        debugfile='debug_cryspy',
    )
    logger = getLogger('cryspy')

    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    else:
        with open('lock_cryspy', 'w'):
            pass    # create vacant file

    # ---------- restart
    if os.path.isfile('cryspy.stat'):
        rin, init_struc_data, rng = cryspy_restart.restart()
    else:
        logger.error('cryspy.stat file does not exist.')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- algo check
    if rin.algo == 'RS':
        plot_RS(rin)
    elif rin.algo == 'EA':
        plot_EA(rin)
    elif rin.algo == 'EA-vc':
        plot_EA_vc(rin)
    else:
        logger.error(f'algo = {rin.algo} is not supported in cryspy-Eplot')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- unlock
    os.remove('lock_cryspy')


# ----------
# ---------- RS
# ----------
def plot_RS(rin):
    '''
    Plot energy for RS.
    '''

    # ---------- logger
    logger = getLogger('cryspy')

    # ---------- load data
    rslt_data = pkl_data.load_rslt()

    # ---------- plot
    fig, ax = plot_energy_RS(rslt_data, rin.ymax, rin.markersize)

    # ---------- save figure
    os.makedirs('./data/energy_plot', exist_ok=True)
    fname = f'./data/energy_plot/energy_RS_{len(rslt_data)}.{rin.fig_format}'
    fig.savefig(fname)
    logger.info(f'Energy plot for RS saved as {fname}')


# ----------
# ---------- EA
# ----------
def plot_EA(rin):
    '''
    Plot energy for EA.
    '''

    # ---------- logger
    logger = getLogger('cryspy')

    # ---------- load data
    rslt_data = pkl_data.load_rslt()

    # ---------- plot
    fig, ax = plot_energy_EA(rslt_data, rin.ymax, rin.markersize)

    # ---------- save figure
    os.makedirs('./data/energy_plot', exist_ok=True)
    gmax = rslt_data['Gen'].max()
    fname = f'./data/energy_plot/energy_EA_gen_{gmax}.{rin.fig_format}'
    fig.savefig(fname)
    logger.info(f'Energy plot for EA saved as {fname}')


# ----------
# ---------- EA-vc
# ----------
def plot_EA_vc(rin):
    '''
    Plot convex hull for EA-vc.
    '''

    # ---------- logger
    logger = getLogger('cryspy')

    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()

    # ---------- generation range
    g_max_avail = max(pd_data.keys())
    # ------ plot_min_gen
    if rin.plot_min_gen is not None:
        if rin.plot_min_gen > g_max_avail:
            logger.error(f'plot_min_gen = {rin.plot_min_gen} is larger than the maximum generation in pd_data')
            logger.error(f'Latest generation in pd_data is {g_max_avail}')
            os.remove('lock_cryspy')
            raise SystemExit(1)
        g_min = rin.plot_min_gen
    else:
        g_min = 1
    # ------ plot_max_gen
    if rin.plot_max_gen is not None:
        if rin.plot_max_gen > g_max_avail:
            logger.error(f'plot_max_gen = {rin.plot_max_gen} is larger than the maximum generation in pd_data')
            logger.error(f'Latest generation in pd_data is {g_max_avail}')
            os.remove('lock_cryspy')
            raise SystemExit(1)
        g_max = rin.plot_max_gen
    else:
        g_max = g_max_avail
    # ------ hull_ref_gen
    if rin.hull_ref_gen is not None:
        if rin.hull_ref_gen > g_max_avail:
            logger.error(f'hull_ref_gen = {rin.hull_ref_gen} is larger than the maximum generation in pd_data')
            logger.error(f'Latest generation in pd_data is {g_max_avail}')
            os.remove('lock_cryspy')
            raise SystemExit(1)
        hull_ref_gen = rin.hull_ref_gen
    elif rin.plot_max_gen is not None:
        hull_ref_gen = rin.plot_max_gen
    else:
        hull_ref_gen = g_max_avail

    # ---------- phase_diagram, hdist
    phase_diagram = pd_data[hull_ref_gen]
    hdist = hdist_data[hull_ref_gen]

    # ---------- filtering generations
    if rin.plot_min_gen is not None or rin.plot_max_gen is not None:
        filtered_rslt = rslt_data[(rslt_data['Gen'] >= g_min) & (rslt_data['Gen'] <= g_max)]
        filtered_ids = filtered_rslt.index.values
    else:
        filtered_ids = None

    # ---------- axis order
    atype = rin.atype
    axis_order = rin.axis_order
    if len(atype) == 2:    # ordering not used for binary system
        if axis_order is None:
            axis_order = 'lr'
    if len(atype) == 3:
        ordering = None
        if axis_order is None:
            axis_order ='tlr'
        ordering = _build_ordering(atype, axis_order)

    # ---------- plot
    if len(rin.atype) == 2:
        fig, ax = draw_convex_hull_binary(
            phase_diagram,
            hdist,
            filtered_ids,
            rin.ymax,
            rin.label_stable,
            rin.vmax,
            rin.bottom_margin,
            rin.markersize,
            axis_order,
        )
    elif len(rin.atype) == 3:
        fig, ax = draw_convex_hull_ternary(
            phase_diagram,
            hdist,
            filtered_ids,
            rin.show_max,
            rin.label_stable,
            rin.vmax,
            rin.markersize,
            ordering,
        )
    else:
        logger.error(f'len(rin.atype) = {len(rin.atype)} is not supported in cryspy-Eplot')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- save figure
    if rin.plot_min_gen is None and rin.plot_max_gen is None:
        # ------ normal mode (latest only)
        fname = f'./data/convex_hull/conv_hull_gen_{g_max_avail}.{rin.fig_format}'
    elif rin.plot_min_gen is None and rin.plot_max_gen is not None:
        # ------ only max specified
        if rin.hull_ref_gen is None:
            # normal mode
            fname = f'./data/convex_hull/conv_hull_gen_{rin.plot_max_gen}.{rin.fig_format}'
        else:
            # with ref gen
            if rin.hull_ref_gen == rin.plot_max_gen:
                # normal mode
                fname = f'./data/convex_hull/conv_hull_gen_{rin.plot_max_gen}.{rin.fig_format}'
            else:
                # different ref gen
                fname = (
                    f'./data/convex_hull/conv_hull_gen_{rin.plot_max_gen}'
                    f'_href_{rin.hull_ref_gen}.{rin.fig_format}'
                )
    else:
        # ------ min exists (range plot)
        href = rin.hull_ref_gen if rin.hull_ref_gen is not None else g_max
        fname = (
            f'./data/convex_hull/conv_hull_gen_{g_min}-{g_max}'
            f'_href_{href}.{rin.fig_format}'
        )
    # ------ save
    os.makedirs('data/convex_hull', exist_ok=True)
    fig.savefig(fname)
    logger.info(f'Convex hull plot saved as {fname}')


def _build_ordering(atype, axis_order):
    """
    Build ordering list for PDPlotter / order_phase_diagram.

    Interpretation of `axis_order`:
        Binary (len(atype) == 2):
            axis_order is a length-2 string consisting of 'l' and 'r'.
            axis_order[i] indicates whether atype[i] is plotted on
            Left or Right of the binary phase diagram.

        Ternary (len(atype) == 3):
            axis_order is a length-3 string consisting of 't', 'l', 'r'.
            axis_order[i] indicates whether atype[i] is plotted on
            Top, Left, or Right of the ternary triangle.

    Example:
        atype = ('Li', 'Ca', 'Cl')
        axis_order = "rtl"
            Li -> Right
            Ca -> Top
            Cl -> Left
        returns ['Ca', 'Cl', 'Li']   # [Up, Left, Right]

    Returns:
        list[str] or None
            For binary: [Left, Right]
            For ternary: [Up, Left, Right]
            For other dimensions: None (ordering not applied)
    """

    # ---------- binary system
    if len(atype) == 2:
        idx_l = axis_order.index('l')
        idx_r = axis_order.index('r')
        left  = _to_special_formula(atype[idx_l])
        right = _to_special_formula(atype[idx_r])
        return [left, right]

    # ---------- ternary system
    elif len(atype) == 3:
        idx_t = axis_order.index('t')
        idx_l = axis_order.index('l')
        idx_r = axis_order.index('r')
        up    = _to_special_formula(atype[idx_t])
        left  = _to_special_formula(atype[idx_l])
        right = _to_special_formula(atype[idx_r])
        return [up, left, right]

    # ---------- higher components
    # ordering control is not supported for >3 components
    return None


def _to_special_formula(sym: str) -> str:
    special_formulas = {
        "O":  "O2",
        "N":  "N2",
        "F":  "F2",
        "Cl": "Cl2",
        "H":  "H2",
    }
    return special_formulas.get(sym, sym)


# ----------
# ---------- main
# ----------
if __name__ == '__main__':
    # ---------- main
    main()