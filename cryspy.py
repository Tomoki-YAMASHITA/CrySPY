#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import CrySPY
from CrySPY.IO import read_input as rin


#---------- initialize
if not os.path.isfile('cryspy.stat'):
    #------ cryspy_init
    stat, init_struc_data, opt_struc_data, rslt_data = CrySPY.start.cryspy_init.initialize()
    if rin.algo == 'RS':
        RS_id_data = CrySPY.start.cryspy_init.RS_init(stat)
    elif rin.algo == 'BO':
        rslt_data, BO_id_data, BO_data= CrySPY.BO.BO_init.initialize(stat, init_struc_data, rslt_data)
    if rin.kpt_flag:
        kpt_data = CrySPY.start.cryspy_init.kpt_init()

#---------- restart
else:
    #------ cpsy_restart
    stat = CrySPY.start.cryspy_restart.restart()

    #------ load data
    init_struc_data = CrySPY.IO.pkl_data.load_init_struc()
    opt_struc_data = CrySPY.IO.pkl_data.load_opt_struc()
    rslt_data = CrySPY.IO.pkl_data.load_rslt()
    if rin.algo == 'RS':
        RS_id_data = CrySPY.IO.pkl_data.load_RS_id()
    elif rin.algo == 'BO':
        BO_id_data = CrySPY.IO.pkl_data.load_BO_id()
        BO_data = CrySPY.IO.pkl_data.load_BO_data()
    if rin.kpt_flag:
        kpt_data = CrySPY.IO.pkl_data.load_kpt()

    #------ append structures
    if len(init_struc_data) < rin.tot_struc:
        init_struc_data = CrySPY.start.cryspy_restart.append_struc(init_struc_data)
    elif rin.tot_struc < len(init_struc_data):
        raise ValueError('tot_struc < len(init_struc_data)')
    #-- BO
    if rin.algo == 'BO':
        if BO_id_data[1] < len(init_struc_data):    # BO_id_data[1] is next_BO_id
            BO_id_data, BO_data = CrySPY.BO.BO_restart.restart(init_struc_data, BO_id_data, BO_data)




#---------- check point 1
if rin.stop_chkpt == 1:
    print('Stop at check point 1')
    raise SystemExit()




#---------- check calc files in ./calc_in
CrySPY.interface.select_code.check_calc_files()




#---------- check point 2
if rin.stop_chkpt == 2:
    print('Stop at check point 2')
    raise SystemExit()




#---------- make working directory
if not os.path.isdir('work{:04d}'.format(rin.njob - 1)):
    for i in range(rin.njob):
        if not os.path.isdir('work{:04d}'.format(i)):
            os.mkdir('work{:04d}'.format(i))

#---------- instantiate Ctrl_job class
jobs = CrySPY.job.ctrl_job.Ctrl_job(stat, init_struc_data, opt_struc_data, rslt_data)
if rin.algo == 'RS':
    jobs.RS_init(RS_id_data)
elif rin.algo == 'BO':
    jobs.BO_init(BO_id_data, BO_data)
if rin.kpt_flag:
    jobs.kpt_init(kpt_data)

#---------- check job status
jobs.check_job()

#---------- control job
print('\n#---------- job status')
for work_id, jstat in enumerate(jobs.job_stat):
    #------ set work_id and work_path
    jobs.work_id = work_id
    jobs.work_path = './work{:04d}/'.format(work_id)

    #------ handle job
    if jstat == 'submitted':
        print('work{:04d}: still queuing or runnning'.format(work_id))
    elif jstat == 'done':
        jobs.handle_done()
    elif jstat == 'skip':
        jobs.ctrl_skip()
    elif jstat == 'else':
        raise ValueError('Wrong job_stat in ' +
                         jobs.work_path +
                         'stat_job. To skip this structure, write "skip" in stat_job line 3')
    elif jstat == 'no_file':
        jobs.ctrl_next_struc()
    else:
        raise ValueError('Unexpected error in '+jobs.work_path+'stat_job')

#---------- BO
if rin.algo == 'BO':
    if jobs.logic_next_gen:
        #------ log and out
        with open('cryspy.out', 'a') as fout:
            fout.write('\nDone generation {}\n\n'.format(jobs.gen))
        print('\nDone generation {}\n'.format(jobs.gen))

        #------ done all structures
        if len(jobs.rslt_data) == rin.tot_struc:
            with open('cryspy.out', 'a') as fout:
                fout.write('\nDone all structures!\n')
            print('\nDone all structures!')
            raise SystemExit()

        #------ check job status
        jobs.check_job()

        #------ next generation
        if 'submitted' not in jobs.job_stat:

            #---------- check point 3
            if rin.stop_chkpt == 3:
                print('Stop at check point 3: BO is ready')
                raise SystemExit()

            jobs.ctrl_next_gen()
