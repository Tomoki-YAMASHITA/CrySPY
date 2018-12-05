#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

from CrySPY.BO import bo_init, bo_restart, bo_next_gen
from CrySPY.EA import ea_init, ea_next_gen, ea_append
from CrySPY.interface import select_code
from CrySPY.job import ctrl_job
from CrySPY.IO import pkl_data
from CrySPY.IO import read_input as rin
from CrySPY.LAQA import laqa_init, laqa_restart, laqa_next_selection
from CrySPY.start import cryspy_init, cryspy_restart


# ---------- lock
if os.path.isfile('lock_cryspy'):
    raise SystemExit('lock_cryspy file exists')
else:
    with open('lock_cryspy', 'w') as f:
        pass    # create vacant file

# ---------- initialize
if not os.path.isfile('cryspy.stat'):
    # ------ cryspy_init
    stat, init_struc_data, rslt_data = cryspy_init.initialize()
    if rin.algo == 'RS':
        cryspy_init.rs_init(stat)
    elif rin.algo == 'BO':
        bo_init.initialize(stat, init_struc_data, rslt_data)
    elif rin.algo == 'LAQA':
        laqa_init.initialize(stat, init_struc_data)
    elif rin.algo == "EA":
        ea_init.initialize(stat, rslt_data)
    if rin.kpt_flag:
        cryspy_init.kpt_init()
    if rin.energy_step_flag:
        cryspy_init.energy_step_init()
    if rin.struc_step_flag:
        cryspy_init.struc_step_init()
    if rin.fs_step_flag:
        cryspy_init.fs_step_init()
    os.remove('lock_cryspy')
    raise SystemExit()

# ---------- restart
else:
    # ------ cpsy_restart
    stat = cryspy_restart.restart()
    # ------ load init_struc_data for appending structures
    init_struc_data = pkl_data.load_init_struc()
    # ------ append structures
    if len(init_struc_data) < rin.tot_struc:
        prev_nstruc = len(init_struc_data)
        init_struc_data = cryspy_restart.append_struc(init_struc_data)
        # -- BO
        if rin.algo == 'BO':
            bo_id_data = pkl_data.load_bo_id()
            bo_data = pkl_data.load_bo_data()
            bo_restart.restart(init_struc_data, bo_id_data, bo_data, prev_nstruc)
        # -- LAQA
        if rin.algo == 'LAQA':
            laqa_id_data = pkl_data.load_laqa_id()
            laqa_data = pkl_data.load_laqa_data()
            laqa_restart.restart(stat, laqa_id_data, laqa_data, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()
    elif rin.tot_struc < len(init_struc_data):
        raise ValueError('tot_struc < len(init_struc_data)')
    # ------ load data
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()
    # ------ append structures by EA
    if rin.append_struc_ea:
        prev_nstruc = len(init_struc_data)
        init_struc_data = ea_append.append_struc(stat, init_struc_data, opt_struc_data, rslt_data)
        # -- BO
        if rin.algo == 'BO':
            bo_id_data = pkl_data.load_bo_id()
            bo_data = pkl_data.load_bo_data()
            bo_restart.restart(init_struc_data, bo_id_data, bo_data, prev_nstruc)
        # -- LAQA
        if rin.algo == 'LAQA':
            laqa_id_data = pkl_data.load_laqa_id()
            laqa_data = pkl_data.load_laqa_data()
            laqa_restart.restart(stat, laqa_id_data, laqa_data, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()
    if rin.algo == 'RS':
        rs_id_data = pkl_data.load_rs_id()
    elif rin.algo == 'BO':
        bo_id_data = pkl_data.load_bo_id()
        bo_data = pkl_data.load_bo_data()
    elif rin.algo == 'LAQA':
        laqa_id_data = pkl_data.load_laqa_id()
        laqa_data = pkl_data.load_laqa_data()
    elif rin.algo == 'EA':
        ea_id_data = pkl_data.load_ea_id()
        # do not have to load ea_data here. ea_data is used only in ea_next_gen.py
    if rin.kpt_flag:
        kpt_data = pkl_data.load_kpt()
    if rin.energy_step_flag:
        energy_step_data = pkl_data.load_energy_step()
    if rin.struc_step_flag:
        struc_step_data = pkl_data.load_struc_step()
    if rin.fs_step_flag:
        fs_step_data = pkl_data.load_fs_step()

#
#
# ---------- check point 1
#
#
if rin.stop_chkpt == 1:
    print('Stop at check point 1')
    os.remove('lock_cryspy')
    raise SystemExit()

# ---------- check calc files in ./calc_in
select_code.check_calc_files()

#
#
# ---------- check point 2
#
#
if rin.stop_chkpt == 2:
    print('Stop at check point 2')
    os.remove('lock_cryspy')
    raise SystemExit()

# ---------- make working directory
if not os.path.isdir('work{:04d}'.format(rin.njob - 1)):
    for i in range(rin.njob):
        if not os.path.isdir('work{:04d}'.format(i)):
            os.mkdir('work{:04d}'.format(i))

# ---------- instantiate Ctrl_job class
jobs = ctrl_job.Ctrl_job(stat, init_struc_data, opt_struc_data, rslt_data)
if rin.algo == 'RS':
    jobs.rs_init(rs_id_data)
elif rin.algo == 'BO':
    jobs.bo_init(bo_id_data, bo_data)
elif rin.algo == 'LAQA':
    jobs.laqa_init(laqa_id_data, laqa_data)
elif rin.algo == 'EA':
    jobs.ea_init(ea_id_data)
if rin.kpt_flag:
    jobs.kpt_data = kpt_data
if rin.energy_step_flag:
    jobs.energy_step_data = energy_step_data
if rin.struc_step_flag:
    jobs.struc_step_data = struc_step_data
if rin.fs_step_flag:
    jobs.fs_step_data = fs_step_data

# ---------- check job status
jobs.check_job()

# ---------- control job
print('\n# ---------- job status')
for work_id, jstat in enumerate(jobs.job_stat):
    # ------ set work_id and work_path
    jobs.work_id = work_id
    jobs.work_path = './work{:04d}/'.format(work_id)

    # ------ set current ID and stage
    jobs.cID = jobs.id_stat[work_id]
    jobs.cstage = jobs.stage_stat[work_id]

    # ------ handle job
    if jstat == 'submitted':
        print('work{:04d}: still queuing or runnning'.format(work_id))
    elif jstat == 'done':
        jobs.handle_done()
    elif jstat == 'skip':
        jobs.ctrl_skip()
        jobs.ctrl_next_struc()
    elif jstat == 'else':
        raise ValueError('Wrong job_stat in ' +
                         jobs.work_path +
                         'stat_job. To skip this structure, write "skip" in stat_job line 3')
    elif jstat == 'no_file':
        jobs.ctrl_next_struc()
    else:
        raise ValueError('Unexpected error in '+jobs.work_path+'stat_job')

# ---------- BO
if rin.algo == 'BO':
    if jobs.logic_next_gen:
        # ------ check job status
        jobs.check_job()
        # ------ next generation
        if set(jobs.job_stat) == {'no_file'}:
            # -- log and out
            with open('cryspy.out', 'a') as fout:
                fout.write('\nDone generation {}\n\n'.format(jobs.gen))
            print('\nDone generation {}\n'.format(jobs.gen))
            # -- done all structures
            if len(jobs.rslt_data) == rin.tot_struc:
                with open('cryspy.out', 'a') as fout:
                    fout.write('\nDone all structures!\n')
                print('\nDone all structures!')
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- check point 3
            if rin.stop_chkpt == 3:
                print('\nStop at check point 3: BO is ready\n')
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- maxgen
            if 0 < rin.maxgen <= jobs.gen:
                print('\nReached maxgen: {}\n'.format(rin.maxgen))
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- BO
            bo_data = (jobs.descriptors, jobs.targets)
            bo_id_data = (jobs.gen, jobs.non_error_id, jobs.id_to_calc, jobs.id_done)
            bo_next_gen.next_gen(jobs.stat, bo_id_data, bo_data)

# ---------- LAQA
elif rin.algo == 'LAQA':
    if jobs.logic_next_selection:
        # ------ check job status
        jobs.check_job()
        # ------ next selection
        if set(jobs.job_stat) == {'no_file'}:
            # -- check point 3
            if rin.stop_chkpt == 3:
                print('\nStop at check point 3: LAQA is ready\n')
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- selection of LAQA
            laqa_id_data = (jobs.id_to_calc, jobs.id_select_hist, jobs.id_done)
            laqa_data = (jobs.tot_step_select, jobs.laqa_step, jobs.laqa_struc,
                         jobs.laqa_energy, jobs.laqa_bias, jobs.laqa_score)
            laqa_next_selection.next_selection(jobs.stat, laqa_id_data, laqa_data)

# ----------  EA
elif rin.algo == 'EA':
    if jobs.logic_next_gen :
        # ------ check job status
        jobs.check_job()
        # ------ next generation:
        if set(jobs.job_stat) == {'no_file'}:
            # -- log and out
            with open('cryspy.out', 'a') as fout:
                fout.write('\nDone generation {}\n\n'.format(jobs.gen))
            print('\nDone generation {}\n'.format(jobs.gen))
            # -- check point 3
            if rin.stop_chkpt == 3:
                print('\nStop at check point 3: EA is ready\n')
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- maxgen
            if 0 < rin.maxgen <= jobs.gen:
                print('\nReached maxgen: {}\n'.format(rin.maxgen))
                os.remove('lock_cryspy')
                raise SystemExit()
            # -- EA
            ea_id_data = (jobs.gen, jobs.next_id, jobs.id_done)
            ea_next_gen.next_gen(jobs.stat, jobs.init_struc_data, jobs.opt_struc_data, jobs.rslt_data, ea_id_data)

# ---------- unlock
os.remove('lock_cryspy')
