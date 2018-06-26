#!/usr/bin/env python
# -*- coding: utf-8 -*-

from VASP import calc_files_vasp, ctrl_job_vasp, collect_vasp
from QE import calc_files_qe, ctrl_job_qe, collect_qe
from soiap import calc_files_soiap, ctrl_job_soiap, collect_soiap
from LAMMPS import calc_files_lammps, ctrl_job_lammps, collect_lammps

from ..IO import read_input as rin


def check_calc_files():
    if rin.calc_code == 'VASP':
        calc_files_vasp.check_input_vasp()
    elif rin.calc_code == 'QE':
        calc_files_qe.check_input_qe()
    elif rin.calc_code == 'soiap':
        calc_files_soiap.check_input_soiap()
    elif rin.calc_code == 'LAMMPS':
        calc_files_lammps.check_input_lammps()
    else:
        raise NotImplementedError('now only VASP, QE, soiap, or LAMMPS')


def next_stage(stage, work_path, *args):
    # args[0] <-- kpt_data
    # args[1] <-- current_ID
    if rin.calc_code == 'VASP':
        skip_flag, kpt_data = ctrl_job_vasp.next_stage_vasp(stage, work_path, args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'QE':
        skip_flag, kpt_data = ctrl_job_qe.next_stage_qe(stage, work_path, args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'soiap':
        skip_flag = ctrl_job_soiap.next_stage_soiap(stage, work_path)
        return skip_flag
    elif rin.calc_code == 'LAMMPS':
        skip_flag = ctrl_job_lammps.next_stage_lammps(stage, work_path)
        return skip_flag
    else:
        raise NotImplementedError('now only VASP, QE, soiap, or LAMMPS')


def collect(current_id, work_path, check_file=None):
    if rin.calc_code == 'VASP':
        opt_struc, energy, magmom, check_opt = \
            collect_vasp.collect_vasp(current_id, work_path, check_file)
    elif rin.calc_code == 'QE':
        opt_struc, energy, magmom, check_opt = \
            collect_qe.collect_qe(current_id, work_path)
    elif rin.calc_code == 'soiap':
        opt_struc, energy, magmom, check_opt = \
            collect_soiap.collect_soiap(current_id, work_path)
    elif rin.calc_code == 'LAMMPS':
        opt_struc, energy, magmom, check_opt = \
            collect_lammps.collect_lammps(current_id, work_path)
    else:
        raise NotImplementedError('now only VASP, QE, soiap, or LAMMPS')

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def next_struc(structure, next_id, work_path, *args):
    # args[0] <-- kpt_data
    if rin.calc_code == 'VASP':
        kpt_data = ctrl_job_vasp.next_struc_vasp(structure, next_id, work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'QE':
        kpt_data = ctrl_job_qe.next_struc_qe(structure, next_id, work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'soiap':
        ctrl_job_soiap.next_struc_soiap(structure, next_id, work_path)
    elif rin.calc_code == 'LAMMPS':
        ctrl_job_lammps.next_struc_lammps(structure, next_id, work_path)
    else:
        raise NotImplementedError('now only VASP, QE, soiap, or LAMMPS')


def clean_calc_files(work_path):
    if rin.calc_code == 'VASP':
        calc_files_vasp.clean_calc_files_vasp(work_path)
    elif rin.calc_code == 'QE':
        calc_files_qe.clean_calc_files_qe(work_path)
    elif rin.calc_code == 'soiap':
        calc_files_soiap.clean_calc_files_soiap(work_path)
    elif rin.calc_code == 'LAMMPS':
        calc_files_lammps.clean_calc_files_lammps(work_path)
    else:
        raise NotImplementedError('now only VASP, QE, soiap, or LAMMPS')


def get_energy_step(energy_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        energy_step_data = collect_vasp.get_energy_step_vasp(energy_step_data, current_id, work_path+'vasprun.xml')
        return energy_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')


def get_struc_step(struc_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        struc_step_data = collect_vasp.get_struc_step_vasp(struc_step_data, current_id, work_path+'vasprun.xml')
        return struc_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')


def get_fs_step(fs_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        fs_step_data = collect_vasp.get_fs_step_vasp(fs_step_data, current_id, work_path+'vasprun.xml')
        return fs_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')
