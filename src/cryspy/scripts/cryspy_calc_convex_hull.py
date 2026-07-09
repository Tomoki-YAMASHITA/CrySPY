#!/usr/bin/env python3
'''
Calculate convex hull script
'''
import argparse
from logging import getLogger
import os

from cryspy.EA.calc_hull import calc_convex_hull
from cryspy.IO import pkl_data
from cryspy.IO.out_results import out_hdist
from cryspy.start import cryspy_restart
from cryspy.util.utility import set_logger


def main():
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('gen', type=int, help='Specify the generation')
    args = parser.parse_args()
    if args.gen < 1:
        raise ValueError('Generation must be greater than or equal to 1')
    gen = args.gen

    # ---------- logger
    set_logger(
        noprint=False,
        debug=False,
        logfile='log_cryspy',
        errfile='err_cryspy',
        debugfile='debug_cryspy',
    )
    logger = getLogger('cryspy')
    logger.info(f'# ---------- cryspy-calc-convex-hull command for generation {gen}')

    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    else:
        with open('lock_cryspy', 'w'):
            pass    # create vacant file

    try:
        # ---------- restart
        if os.path.isfile('cryspy.stat'):
            rin, init_struc_data, rng = cryspy_restart.restart()
        else:
            logger.error('cryspy.stat file does not exist.')
            raise SystemExit(1)

        # ---------- algo check
        if rin.algo not in ['EA-vc']:
            logger.error(
                f'algo = {rin.algo} is not supported in '
                'cryspy-calc-convex-hull'
            )
            raise SystemExit(1)

        # ---------- load data
        logger.info('# ------ Load data')
        rslt_data = pkl_data.load_rslt()
        nat_data = pkl_data.load_nat_data()
        pd_data = pkl_data.load_pd_data()
        hdist_data = pkl_data.load_hdist_data()
        ea_info = pkl_data.load_ea_info()

        # ---------- gen validity check
        g_max_avail = int(rslt_data['Gen'].max())
        if gen > g_max_avail:
            logger.error(
                f'gen = {gen} is larger than the latest available generation '
                f'(latest = {g_max_avail})'
            )
            raise SystemExit(1)

        # ---------- min_comp/max_comp
        min_comp = ea_info.loc[ea_info['Gen'] == gen, 'min_comp'].iloc[-1]
        max_comp = ea_info.loc[ea_info['Gen'] == gen, 'max_comp'].iloc[-1]

        # ---------- calc convex hull
        logger.info('# ------ Calculate convex hull')
        phase_diagram, hdist = calc_convex_hull(
            atype=rin.atype,
            gen=gen,
            ref_energies=rin.ref_energies,
            rslt_data=rslt_data,
            nat_data=nat_data,
            ymax=rin.ymax,
            show_max=rin.show_max,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            bottom_margin=rin.bottom_margin,
            markersize=rin.markersize,
            fig_format=rin.fig_format,
            emax_ea=rin.emax_ea,
            emin_ea=rin.emin_ea,
            axis_order=rin.axis_order,
            min_comp=min_comp,
            max_comp=max_comp,
            show_comp_window=rin.show_comp_window,
        )
        logger.info('# ------ Save data')
        out_hdist(gen, hdist, nat_data)
        hdist_data[gen] = hdist
        pd_data[gen] = phase_diagram
        pkl_data.save_hdist_data(hdist_data)
        pkl_data.save_pd_data(pd_data)

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()