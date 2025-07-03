#!/usr/bin/env python3
'''
Energy plot script
'''

from logging import getLogger
import os

from cryspy.IO import pkl_data
from cryspy.start import cryspy_restart
from cryspy.util.utility import set_logger
from cryspy.util.visual_util import draw_convex_hull_binary, draw_convex_hull_ternary


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
        rin, init_struc_data = cryspy_restart.restart()
    else:
        logger.error('cryspy.stat file does not exist.')
        os.remove('lock_cryspy')

    # ---------- algo check
    if rin.algo not in ['EA-vc']:
        logger.error(f'algo = {rin.algo} is not supported in cryspy-Eplot')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()

    # ---------- current generation
    if rin.cgen is None:
        cgen = max(pd_data.keys())
    else:
        cgen = rin.cgen
        if cgen > max(pd_data.keys()):
            logger.error(f'cgen = {cgen} is larger than the maximum generation in pd_data')
            logger.error(f'Latest generation in pd_data is {max(pd_data.keys())}')
            os.remove('lock_cryspy')
            raise SystemExit(1)
    c_rslt = rslt_data[rslt_data['Gen'] == cgen]
    cgen_ids = c_rslt.index.values    # current IDs [array]
    pd = pd_data[cgen]
    hdist = hdist_data[cgen]

    # ---------- plot
    if len(rin.atype) == 2:
        fig, ax = draw_convex_hull_binary(pd, hdist, cgen_ids, rin.show_max, rin.label_stable, rin.vmax, rin.bottom_margin)
    elif len(rin.atype) == 3:
        fig, ax = draw_convex_hull_ternary(pd, hdist, cgen_ids, rin.show_max, rin.label_stable, rin.vmax)
    else:
        logger.error(f'len(rin.atype) = {len(rin.atype)} is not supported in cryspy-Eplot')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- save figure
    fig.savefig(f'./data/convex_hull/conv_hull_gen_{cgen}.{rin.fig_format}', bbox_inches='tight')
    logger.info(f'Convex hull plot saved as ./data/convex_hull/conv_hull_gen_{cgen}.{rin.fig_format}')

    # ---------- unlock
    os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()