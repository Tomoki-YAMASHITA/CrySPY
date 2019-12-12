#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

from CrySPY.BO import bo_init, bo_restart
from CrySPY.EA import ea_init, ea_append
from CrySPY.interface import select_code
from CrySPY.job import ctrl_job
from CrySPY.IO import pkl_data
from CrySPY.IO import read_input as rin
from CrySPY.LAQA import laqa_init, laqa_restart
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
            bo_data = pkl_data.load_bo_data()
            bo_restart.restart(init_struc_data, bo_data, prev_nstruc)
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

# ---------- handle job
jobs.handle_job()

# ---------- unlock
os.remove('lock_cryspy')
