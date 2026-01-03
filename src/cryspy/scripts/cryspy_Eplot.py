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
    get_generation_range,
    build_ordering,
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
    try:
        g_min, g_max, g_ref = get_generation_range(
            plot_min_gen=rin.plot_min_gen,
            plot_max_gen=rin.plot_max_gen,
            hull_ref_gen=rin.hull_ref_gen,
            g_max_avail=g_max_avail,
        )
    except ValueError as e:
        logger.error(str(e))
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- phase_diagram, hdist
    phase_diagram = pd_data[g_ref]
    hdist = hdist_data[g_ref]

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
        if axis_order is None:
            axis_order ='tlr'
        ordering = build_ordering(atype, axis_order)

    # ---------- plot
    if len(rin.atype) == 2:
        fig, ax = draw_convex_hull_binary(
            phase_diagram=phase_diagram,
            hdist=hdist,
            filtered_ids=filtered_ids,
            ymax=rin.ymax,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            bottom_margin=rin.bottom_margin,
            markersize=rin.markersize,
            axis_order=axis_order,
        )
    elif len(rin.atype) == 3:
        fig, ax = draw_convex_hull_ternary(
            phase_diagram=phase_diagram,
            hdist=hdist,
            filtered_ids=filtered_ids,
            show_max=rin.show_max,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            markersize=rin.markersize,
            ordering=ordering,
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


# ----------
# ---------- main
# ----------
if __name__ == '__main__':
    # ---------- main
    main()