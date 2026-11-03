#!/usr/bin/env python3
'''
Main script
'''

import argparse
from logging import getLogger
import os

from cryspy import __version__
from cryspy.start import cryspy_init, cryspy_restart
from cryspy.util.utility import (
    backup_cryspy,
    clean_cryspy,
    set_logger,
)

# ---------- import later
# from cryspy.job.ctrl_job import Ctrl_job
# from cryspy.interface import select_code

def main():
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--backup', help='backup data', action='store_true')
    parser.add_argument('-c', '--clean', help='clean data', action='store_true')
    parser.add_argument('-g', '--debug', help='debug', action='store_true')
    parser.add_argument('-n', '--noprint', help='not printing to the console', action='store_true')
    args = parser.parse_args()

    # ---------- logger
    set_logger(
        noprint=args.noprint,
        debug=args.debug,
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
        # ---------- check mode
        if os.path.isfile('data/db_data/rslt_data.db'):
            logger.error(
                'data/db_data/rslt_data.db exists for '
                'high-throughput mode'
            )
            raise SystemExit(1)

        # ---------- banner
        logger.info(f'\n\n\nCrySPY {__version__}\n\n')

        # ---------- backup option
        if args.backup:
            backup_cryspy(manual=True)
            raise SystemExit()

        # ---------- clean option
        if args.clean:
            clean_cryspy()
            raise SystemExit()

        # ---------- initialize
        if not os.path.isfile('cryspy.stat'):
            cryspy_init.initialize()
            raise SystemExit()

        # ---------- restart
        else:
            rin, init_struc_data, rng = cryspy_restart.restart()

        # ---------- check calc files in ./calc_in
        from cryspy.interface import select_code
        select_code.check_calc_files(rin)

        # ---------- check point 1
        if rin.stop_chkpt == 1:
            logger.info('Stop at check point 1')
            raise SystemExit()

        # ---------- mkdir work/fin
        os.makedirs('work/fin', exist_ok=True)

        # ---------- instantiate Ctrl_job class
        from cryspy.job.ctrl_job import Ctrl_job
        jobs = Ctrl_job(rin, init_struc_data, rng)

        # ---------- check job status
        jobs.check_job()

        # ---------- handle job
        jobs.handle_job()

        # ---------- recheck for skip and done
        if jobs.id_queueing:
            cnt_recheck = 0
            while jobs.recheck:
                cnt_recheck += 1
                jobs.recheck = False    # True --> False
                logger.info(f'\n\nrecheck {cnt_recheck}\n')
                jobs.check_job()
                jobs.handle_job()

        # ---------- next selection or generation
        if not (jobs.id_queueing or jobs.id_running):
            # ---------- next selection or generation
            if rin.algo in ['BO', 'LAQA', 'EA', 'EA-vc']:
                jobs.next_sg(args.noprint)
            # ---------- for RS
            else:
                logger.info('\nDone all structures!')

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()
