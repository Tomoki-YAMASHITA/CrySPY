#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import VASP
from ..IO import read_input as rin


def check_calc_files():
    if rin.calc_code == 'VASP':
        VASP.calc_files_vasp.check_input_vasp()
    elif rin.calc_code == 'QE':
        raise SystemExit('under construction')
    else:
        raise SystemExit('now only VASP')


def next_stage(stage, work_path):
    if rin.calc_code == 'VASP':
        VASP.ctrl_job_vasp.next_stage_vasp(stage, work_path)
    elif rin.calc_code == 'QE':
        raise SystemExit('under construction')
    else:
        raise SystemExit('now only VASP')


def collect(current_id, work_path):
    if rin.calc_code == 'VASP':
        opt_struc, energy, magmom, check_opt = \
            VASP.collect_vasp.collect_vasp(current_id, work_path)
    elif rin.calc_code == 'QE':
        raise SystemExit('under construction')
    else:
        raise SystemExit('now only VASP')

    #---------- return
    return opt_struc, energy, magmom, check_opt


def next_struc(init_struc_data, next_id, work_path):
    if rin.calc_code == 'VASP':
        VASP.ctrl_job_vasp.next_struc_vasp(init_struc_data, next_id, work_path)
    elif rin.calc_code == 'QE':
        raise SystemExit('under construction')
    else:
        raise SystemExit('now only VASP')


def clean_calc_files(work_path):
    if rin.calc_code == 'VASP':
        VASP.calc_files_vasp.clean_calc_files_vasp(work_path)
    elif rin.calc_code == 'QE':
        raise SystemExit('under construction')
    else:
        raise SystemExit('now only VASP')
