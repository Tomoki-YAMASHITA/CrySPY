#!/usr/bin/env python3
'''
High-throughput mode.
'''
import argparse
from logging import getLogger
import os

from cryspy import __version__
from cryspy.high_throughput.db.ea import select_latest_ea_generation
from cryspy.high_throughput.db.sqlite import connect_db
from cryspy.high_throughput.start.check_input import check_calc_file
from cryspy.high_throughput.start.initialize import initialize
from cryspy.high_throughput.start.restart import restart
from cryspy.high_throughput.worker.controller_ea import (
    generate_next_generation_ea,
)
from cryspy.high_throughput.worker.controller_opt import (
    RunStatus,
    launch_workers_opt,
)
from cryspy.util.utility import (
    backup_cryspy,
    clean_cryspy,
    set_logger,
)


def main():
    # ---------- argparse
    parser = argparse.ArgumentParser(description='CrySPY high-throughput mode')
    parser.add_argument('-b', '--backup', help='backup data', action='store_true')
    parser.add_argument('-c', '--clean', help='clean data', action='store_true')
    parser.add_argument('-g', '--debug', help='debug', action='store_true')
    parser.add_argument('-p', '--print', help='print logs to the console', action='store_true')
    args = parser.parse_args()

    # ---------- logger
    set_logger(
        noprint=not args.print,
        debug=args.debug,
        logfile='log_cryspy-ht',
        errfile='err_cryspy-ht',
        debugfile='debug_cryspy-ht',
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
        if os.path.isfile('cryspy.stat'):
            logger.error('cryspy.stat exists for normal mode')
            raise SystemExit(1)

        logger.info(
            f'\n\n\nCrySPY high-throughput mode '
            f'{__version__}\n\n'
        )

        # ---------- backup option
        if args.backup:
            backup_cryspy(ht=True, manual=True)
            raise SystemExit()

        # ---------- clean option
        if args.clean:
            clean_cryspy(ht=True)
            raise SystemExit()

        # ---------- initialize or restart
        if not os.path.isfile('data/db_data/rslt_data.db'):
            rin = initialize()
        else:
            rin = restart()

        # ---------- check point 2
        if rin.stop_chkpt == 2:
            logger.info('Stop at check point 2')
            raise SystemExit()

        # ---------- check calculation file
        check_calc_file(rin)

        # ---------- main execution loop
        while True:
            # ------ structure optimization with multiple workers
            run_status = launch_workers_opt(rin)

            # ------ stop
            if run_status == RunStatus.STOPPED:
                break

            # ------ random search
            if rin.algo == 'RS':
                logger.info('\nDone all structures!')
                break

            # ------ current generation
            with connect_db() as conn:
                generation = select_latest_ea_generation(conn)
            if generation is None:
                logger.error('EA generation information was not found')
                raise SystemExit(1)
            logger.info(f'Done generation {generation}')

            # ------ check point 3
            if rin.stop_chkpt == 3:
                logger.info(
                    '\nStop at check point 3: EA is ready'
                )
                break

            # ------ maxgen_ea
            if 0 < rin.maxgen_ea <= generation:
                logger.info(
                    f'\nReached maxgen_ea: {rin.maxgen_ea}'
                )
                break

            # ------ backup
            if 0 < rin.backup_interval and generation % rin.backup_interval == 0:
                backup_cryspy(ht=True)

            # ------ next generation
            generate_next_generation_ea(
                rin,
                generation,
            )

    finally:
        # ---------- unlock
        os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()
