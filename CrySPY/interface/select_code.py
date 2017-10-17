#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import VASP
from . import QE
from . import soiap
from ..IO import read_input as rin


def check_calc_files():
    if rin.calc_code == 'VASP':
        VASP.calc_files_vasp.check_input_vasp()
    elif rin.calc_code == 'QE':
        QE.calc_files_qe.check_input_qe()
    elif rin.calc_code == 'soiap':
        soiap.calc_files_soiap.check_input_soiap()
    else:
        raise SystemExit('now only VASP, QE, or soiap')


def next_stage(stage, work_path, *args):
    # args[0] <-- kpt_data
    # args[1] <-- current_ID
    if rin.calc_code == 'VASP':
        skip_flag, kpt_data = VASP.ctrl_job_vasp.next_stage_vasp(stage, work_path, args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'QE':
        skip_flag, kpt_data = QE.ctrl_job_qe.next_stage_qe(stage, work_path, args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'soiap':
        skip_flag = soiap.ctrl_job_soiap.next_stage_soiap(stage, work_path)
        return skip_flag
    else:
        raise SystemExit('now only VASP, QE, or soiap')


def collect(current_id, work_path):
    if rin.calc_code == 'VASP':
        opt_struc, energy, magmom, check_opt = \
            VASP.collect_vasp.collect_vasp(current_id, work_path)
    elif rin.calc_code == 'QE':
        opt_struc, energy, magmom, check_opt = \
            QE.collect_qe.collect_qe(current_id, work_path)
    elif rin.calc_code == 'soiap':
        opt_struc, energy, magmom, check_opt = \
            soiap.collect_soiap.collect_soiap(current_id, work_path)
    else:
        raise SystemExit('now only VASP, QE, or soiap')

    #---------- return
    return opt_struc, energy, magmom, check_opt


def next_struc(init_struc_data, next_id, work_path, *args):
    # args[0] <-- kpt_data
    if rin.calc_code == 'VASP':
        kpt_data = VASP.ctrl_job_vasp.next_struc_vasp(init_struc_data, next_id, work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'QE':
        kpt_data = QE.ctrl_job_qe.next_struc_qe(init_struc_data, next_id, work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'soiap':
        soiap.ctrl_job_soiap.next_struc_soiap(init_struc_data, next_id, work_path)
    else:
        raise SystemExit('now only VASP, QE, or soiap')


def clean_calc_files(work_path):
    if rin.calc_code == 'VASP':
        VASP.calc_files_vasp.clean_calc_files_vasp(work_path)
    elif rin.calc_code == 'QE':
        QE.calc_files_qe.clean_calc_files_qe(work_path)
    elif rin.calc_code == 'soiap':
        soiap.calc_files_soiap.clean_calc_files_soiap(work_path)
    else:
        raise SystemExit('now only VASP, QE, or soiap')
