#!/usr/bin/env python3
'''
Main script
'''

import os

from CrySPY.interface import select_code
from CrySPY.job.ctrl_job import Ctrl_job
from CrySPY.IO import read_input as rin
from CrySPY.start import cryspy_init, cryspy_restart


def main():
    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        raise SystemExit('lock_cryspy file exists')
    else:
        with open('lock_cryspy', 'w') as f:
            pass    # create vacant file

    # ---------- initialize
    if not os.path.isfile('cryspy.stat'):
        cryspy_init.initialize()
        os.remove('lock_cryspy')
        raise SystemExit()
    # ---------- restart
    else:
        stat, init_struc_data = cryspy_restart.restart()

    # ---------- check point 1
    if rin.stop_chkpt == 1:
        print('Stop at check point 1')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- check calc files in ./calc_in
    select_code.check_calc_files()

    # ---------- mkdir work/fin
    os.makedirs('work/fin', exist_ok=True)

    # ---------- instantiate Ctrl_job class
    jobs = Ctrl_job(stat, init_struc_data)

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
            print('\n\n recheck {}\n'.format(cnt_recheck))
            jobs.check_job()
            jobs.handle_job()

    if not (jobs.id_queueing or jobs.id_running):
        # ---------- next selection or generation
        if rin.algo in ['BO', 'LAQA', 'EA']:
            jobs.next_sg()
        # ---------- for RS
        else:
            with open('cryspy.out', 'a') as fout:
                fout.write('\nDone all structures!\n')
                print('Done all structures!')

    # ---------- unlock
    os.remove('lock_cryspy')


if __name__ == '__main__':
    main()
