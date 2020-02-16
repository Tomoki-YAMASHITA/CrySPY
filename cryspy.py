#!/usr/bin/env python
'''
Main script
'''

import os

from CrySPY.interface import select_code
from CrySPY.job import ctrl_job
from CrySPY.IO import read_input as rin
from CrySPY.start import cryspy_init, cryspy_restart


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

# os.makedirs('work/fin', exist_ok=True)    # python3.2 or later
if not os.path.isdir('work/fin'):
    os.makedirs('work/fin')

# ---------- instantiate Ctrl_job class
jobs = ctrl_job.Ctrl_job(stat, init_struc_data)

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

# ---------- next selection or generation
if rin.algo in ['BO', 'LAQA', 'EA']:
    if not (jobs.id_queueing or jobs.id_running):
        jobs.next_sg()

# ---------- unlock
os.remove('lock_cryspy')
