#!/usr/bin/env python3
'''
Energy plot script
'''

import argparse
from logging import getLogger
import os

from cryspy.IO import pkl_data, read_input
from cryspy.util.utility import set_logger
from cryspy.util.visual_util import (
    plot_energy_RS,
    plot_energy_EA,
    draw_convex_hull_binary,
    draw_convex_hull_ternary,
    prepare_EA_data,
    prepare_hull_data,
)


def main():

    # ---------- args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c',
        '--composition-window',
        action='store_true',
        help='Plot composition window from cryspy.in.',
    )
    args = parser.parse_args()

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

    try:
        # ---------- read input
        if os.path.isfile('cryspy.in'):
            rin = read_input.ReadInput('cryspy.in')
        else:
            logger.error('cryspy.in file does not exist.')
            raise SystemExit(1)

        # ---------- composition window
        if args.composition_window:
            plot_composition_window_input(rin)
            return

        # ---------- algo check
        if rin.algo == 'RS':
            plot_RS(rin)
        elif rin.algo == 'EA':
            plot_EA(rin)
        elif rin.algo == 'EA-vc':
            plot_EA_vc(rin)
        else:
            logger.error(
                f'algo = {rin.algo} is not supported in cryspy-Eplot'
            )
            raise SystemExit(1)

    finally:
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

    # ---------- prepare data
    try:
        (
            filtered_rslt,
            max_indx,
            g_ref,
            g_min,
            g_max,
        ) = prepare_EA_data(
            rslt_data=None,
            ref_gen=rin.ref_gen,
            plot_min_gen=rin.plot_min_gen,
            plot_max_gen=rin.plot_max_gen,
        )
    except ValueError as e:
        logger.error(str(e))
        raise SystemExit(1)

    # ---------- plot
    fig, ax = plot_energy_EA(filtered_rslt, rin.ymax, rin.markersize, max_indx=max_indx)

    # ---------- save figure (use effective generation values)
    if g_min == 1 and g_max == g_ref:
        # ------ normal plot as of generation g_ref
        fname = f'./data/energy_plot/energy_EA_gen_{g_ref}.{rin.fig_format}'
    else:
        # ------ range control
        fname = (
            f'./data/energy_plot/energy_EA_gen_{g_min}-{g_max}'
            f'_ref_{g_ref}.{rin.fig_format}'
        )
    # ------ save figure
    os.makedirs('./data/energy_plot', exist_ok=True)
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

    # ---------- prepare data
    try:
        (
            phase_diagram,
            hdist,
            min_comp,
            max_comp,
            filtered_ids,
            g_ref,
            g_min,
            g_max,
        ) = prepare_hull_data(
                rslt_data=None,
                ea_info=None,
                pd_data=None,
                hdist_data=None,
                ref_gen=rin.ref_gen,
                plot_min_gen=rin.plot_min_gen,
                plot_max_gen=rin.plot_max_gen,
                ref_gen_comp=rin.ref_gen_comp,
            )
    except ValueError as e:
        logger.error(str(e))
        raise SystemExit(1)

    # ---------- axis order
    atype = rin.atype
    axis_order = rin.axis_order
    if len(atype) == 2:
        if axis_order is None:
            axis_order = 'lr'
    if len(atype) == 3:
        if axis_order is None:
            axis_order ='tlr'

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
            min_comp=min_comp,
            max_comp=max_comp,
            show_comp_window=rin.show_comp_window,
        )
    elif len(rin.atype) == 3:
        fig, ax = draw_convex_hull_ternary(
            atype=atype,
            phase_diagram=phase_diagram,
            hdist=hdist,
            filtered_ids=filtered_ids,
            show_max=rin.show_max,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            markersize=rin.markersize,
            axis_order=axis_order,
            min_comp=min_comp,
            max_comp=max_comp,
            show_comp_window=rin.show_comp_window,
        )
    else:
        logger.error(f'len(rin.atype) = {len(rin.atype)} is not supported in cryspy-Eplot')
        raise SystemExit(1)

    # ---------- save figure (use effective generation values)
    if g_min == 1 and g_max == g_ref:
        # ------ normal plot as of generation g_ref
        fname = f'./data/convex_hull/conv_hull_gen_{g_ref}.{rin.fig_format}'
    else:
        # ------ range control
        fname = (
            f'./data/convex_hull/conv_hull_gen_{g_min}-{g_max}'
            f'_ref_{g_ref}.{rin.fig_format}'
        )
    # ------ save
    os.makedirs('data/convex_hull', exist_ok=True)
    fig.savefig(fname)
    logger.info(f'Convex hull plot saved as {fname}')


#
# ----------
# ---------- composition window from cryspy.in
# ----------
def plot_composition_window_input(rin):
    '''
    Plot composition window from cryspy.in.
    '''

    # ---------- logger
    logger = getLogger('cryspy')

    # ---------- check
    if rin.algo != 'EA-vc':
        logger.error('Composition window plot is available only for EA-vc.')
        raise SystemExit(1)

    if rin.min_comp is None and rin.max_comp is None:
        logger.error('min_comp and max_comp are not set in cryspy.in.')
        raise SystemExit(1)

    if len(rin.atype) not in (2, 3):
        logger.error(
            f'len(rin.atype) = {len(rin.atype)} is not supported '
            'for composition window plot'
        )
        raise SystemExit(1)

    # ---------- plot
    from cryspy.util.visual_util import save_composition_window
    save_composition_window(
        atype=rin.atype,
        gen='input',
        min_comp=rin.min_comp,
        max_comp=rin.max_comp,
        fig_format=rin.fig_format,
        axis_order=rin.axis_order,
    )

    # ---------- log
    logger.info(
        f'Composition window saved as '
        f'./data/convex_hull/composition_window_input.{rin.fig_format}'
    )


# ----------
# ---------- main
# ----------
if __name__ == '__main__':
    # ---------- main
    main()