#!/usr/bin/env python3
'''
High-throughput mode.
'''
import argparse
from logging import getLogger
import os

from cryspy.high_throughput.start.initialize import initialize
from cryspy.high_throughput.start.restart import restart
from cryspy.high_throughput.worker.run import launch_workers
from cryspy.util.utility import set_logger


def main():
    # ---------- argparse
    parser = argparse.ArgumentParser(description='CrySPY high-throughput mode')
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

    # ---------- initialize or restart
    if not os.path.isfile('data/db_data/rslt_data.db'):
        rin = initialize()
    else:
        rin = restart()

    # ---------- structure optimization with multiple workers
    launch_workers(rin)


if __name__ == '__main__':
    # ---------- main
    main()