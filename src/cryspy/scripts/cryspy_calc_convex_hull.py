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

    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    else:
        with open('lock_cryspy', 'w'):
            pass    # create vacant file

    # ---------- restart
    if os.path.isfile('cryspy.stat'):
        rin, init_struc_data = cryspy_restart.restart()
    else:
        logger.error('cryspy.stat file does not exist.')
        os.remove('lock_cryspy')

    # ---------- algo check
    if rin.algo not in ['EA-vc']:
        logger.error(f'algo = {rin.algo} is not supported in cryspy-calc-convex-hull')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- load data
    logger.info('# ---------- Load data')
    rslt_data = pkl_data.load_rslt()
    nat_data = pkl_data.load_nat_data()
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()

    # ---------- calc convex hull
    logger.info('# ---------- Calculate convex hull')
    pd, hdist = calc_convex_hull(
        atype=rin.atype,
        gen=gen,
        end_point=rin.end_point,
        rslt_data=rslt_data,
        nat_data=nat_data,
        show_max=rin.show_max,
        label_stable=rin.label_stable,
        vmax=rin.vmax,
        bottom_margin=rin.bottom_margin,
        fig_format=rin.fig_format,
        emax_ea=rin.emax_ea,
        emin_ea=rin.emin_ea,
    )
    logger.info('# ---------- Save data')
    out_hdist(gen, hdist, nat_data)
    hdist_data[gen] = hdist
    pd_data[gen] = pd
    pkl_data.save_hdist_data(hdist_data)
    pkl_data.save_pd_data(pd_data)

    # ---------- unlock
    os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()