#!/usr/bin/env python3
'''
Main script
'''

import argparse
from logging import getLogger
import os

from cryspy.start import cryspy_init, cryspy_restart
from cryspy.util.utility import set_logger, backup_cryspy, clean_cryspy

# ---------- import later
# from mpi4py import MPI
# from cryspy.job.ctrl_job import Ctrl_job
# from cryspy.interface import select_code

def main():
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--backup', help='backup data', action='store_true')
    parser.add_argument('-c', '--clean', help='clean data', action='store_true')
    parser.add_argument('-g', '--debug', help='debug', action='store_true')
    parser.add_argument('-n', '--noprint', help='not printing to the console', action='store_true')
    parser.add_argument('-p', '--parallel', help='with MPI', action='store_true')
    args = parser.parse_args()

    # ########## MPI start
    # ---------- MPI
    if args.parallel:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        mpi_rank = comm.Get_rank()
        mpi_size = comm.Get_size()
    else:
        # ------ normal run
        comm = None
        mpi_rank = 0
        mpi_size = 1

    # ---------- logger
    set_logger(
        noprint=args.noprint,
        debug=args.debug,
        logfile='log_cryspy',
        errfile='err_cryspy',
        debugfile='debug_cryspy',
    )
    logger = getLogger('cryspy')

    # ---------- backup option
    if args.backup:
        if mpi_rank == 0:
            backup_cryspy()
        raise SystemExit()

    # ---------- clean option
    if args.clean:
        if mpi_rank == 0:
            clean_cryspy()
        raise SystemExit()

    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        if mpi_rank == 0:
            logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    else:
        if mpi_size > 1:
            comm.barrier()
        if mpi_rank == 0:
            with open('lock_cryspy', 'w'):
                pass    # create vacant file
        else:
            pass

    # ---------- initialize
    if not os.path.isfile('cryspy.stat'):
        cryspy_init.initialize(comm, mpi_rank, mpi_size)
        if mpi_rank == 0:
            os.remove('lock_cryspy')
        raise SystemExit()
    # ---------- restart
    else:
        # only stat and init_struc_data in rank0 are important
        rin, init_struc_data = cryspy_restart.restart(comm, mpi_rank, mpi_size)
    # ########## MPI end

    if mpi_rank == 0:
        # ---------- check point 1
        if rin.stop_chkpt == 1:
            logger.info('Stop at check point 1')
            os.remove('lock_cryspy')
            raise SystemExit()

        # ---------- check calc files in ./calc_in
        from cryspy.interface import select_code
        select_code.check_calc_files(rin)

        # ---------- mkdir work/fin
        os.makedirs('work/fin', exist_ok=True)

        # ---------- instantiate Ctrl_job class
        from cryspy.job.ctrl_job import Ctrl_job
        jobs = Ctrl_job(rin, init_struc_data)

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

        # ---------- unlock
        os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()
