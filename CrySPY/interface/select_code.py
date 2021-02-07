'''
Select an optimizer
    - VASP
    - QE
    - OMX
    - soiap
    - LAMMPS
'''

from .VASP import calc_files_vasp, ctrl_job_vasp, collect_vasp
from .QE import calc_files_qe, ctrl_job_qe, collect_qe
from .soiap import calc_files_soiap, ctrl_job_soiap, collect_soiap
from .LAMMPS import calc_files_lammps, ctrl_job_lammps, collect_lammps
from .OMX import calc_files_OMX, ctrl_job_OMX, collect_OMX

from ..IO import read_input as rin


def check_calc_files():
    if rin.calc_code == 'VASP':
        calc_files_vasp.check_input_vasp()
    elif rin.calc_code == 'QE':
        calc_files_qe.check_input_qe()
    elif rin.calc_code == 'OMX':
        calc_files_OMX.check_input_OMX()
    elif rin.calc_code == 'soiap':
        calc_files_soiap.check_input_soiap()
    elif rin.calc_code == 'LAMMPS':
        calc_files_lammps.check_input_lammps()
    else:
        raise NotImplementedError('now only VASP, QE, OMX, soiap, or LAMMPS')


def next_stage(stage, work_path, *args):
    # args[0] <-- kpt_data
    # args[1] <-- current_id
    if rin.calc_code == 'VASP':
        skip_flag, kpt_data = ctrl_job_vasp.next_stage_vasp(stage, work_path,
                                                            args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'QE':
        skip_flag, kpt_data = ctrl_job_qe.next_stage_qe(stage, work_path,
                                                        args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'OMX':
        skip_flag, kpt_data = ctrl_job_OMX.next_stage_OMX(stage, work_path,
                                                          args[0], args[1])
        return skip_flag, kpt_data
    elif rin.calc_code == 'soiap':
        skip_flag = ctrl_job_soiap.next_stage_soiap(stage, work_path)
        return skip_flag
    elif rin.calc_code == 'LAMMPS':
        skip_flag = ctrl_job_lammps.next_stage_lammps(stage, work_path)
        return skip_flag
    else:
        raise NotImplementedError('now only VASP, QE, OMX, soiap, or LAMMPS')


def collect(current_id, work_path):
    if rin.calc_code == 'VASP':
        opt_struc, energy, magmom, check_opt = \
            collect_vasp.collect_vasp(current_id, work_path)
    elif rin.calc_code == 'QE':
        opt_struc, energy, magmom, check_opt = \
            collect_qe.collect_qe(current_id, work_path)
    elif rin.calc_code == 'OMX':
        opt_struc, energy, magmom, check_opt = \
            collect_OMX.collect_OMX(current_id, work_path)
    elif rin.calc_code == 'soiap':
        opt_struc, energy, magmom, check_opt = \
            collect_soiap.collect_soiap(current_id, work_path)
    elif rin.calc_code == 'LAMMPS':
        opt_struc, energy, magmom, check_opt = \
            collect_lammps.collect_lammps(current_id, work_path)
    else:
        raise NotImplementedError('now only VASP, QE, OMX, soiap, or LAMMPS')

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def next_struc(structure, current_id, work_path, *args):
    # args[0] <-- kpt_data
    if rin.calc_code == 'VASP':
        kpt_data = ctrl_job_vasp.next_struc_vasp(structure, current_id,
                                                 work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'QE':
        kpt_data = ctrl_job_qe.next_struc_qe(structure, current_id,
                                             work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'OMX':
        kpt_data = ctrl_job_OMX.next_struc_OMX(structure, current_id,
                                               work_path, args[0])
        return kpt_data
    elif rin.calc_code == 'soiap':
        ctrl_job_soiap.next_struc_soiap(structure, current_id, work_path)
    elif rin.calc_code == 'LAMMPS':
        ctrl_job_lammps.next_struc_lammps(structure, current_id, work_path)
    else:
        raise NotImplementedError('now only VASP, QE, OMX, soiap, or LAMMPS')


def get_energy_step(energy_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        energy_step_data = collect_vasp.get_energy_step_vasp(
            energy_step_data, current_id, work_path+'vasprun.xml')
        return energy_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        energy_step_data = collect_soiap.get_energy_step_soiap(
            energy_step_data, current_id, work_path+'log.tote')
        return energy_step_data
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')


def get_struc_step(struc_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        struc_step_data = collect_vasp.get_struc_step_vasp(
            struc_step_data, current_id, work_path+'vasprun.xml')
        return struc_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        struc_step_data = collect_soiap.get_struc_step_soiap(
            struc_step_data, current_id, work_path+'log.struc')
        return struc_step_data
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')


def get_force_step(force_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        force_step_data = collect_vasp.get_force_step_vasp(
            force_step_data, current_id, work_path+'vasprun.xml')
        return force_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        force_step_data = collect_soiap.get_force_step_soiap(
            force_step_data, current_id, work_path+'log.frc')
        return force_step_data
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')


def get_stress_step(stress_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        stress_step_data = collect_vasp.get_stress_step_vasp(
            stress_step_data, current_id, work_path+'vasprun.xml')
        return stress_step_data
    elif rin.calc_code == 'QE':
        raise NotImplementedError('now only VASP')
    elif rin.calc_code == 'soiap':
        stress_step_data = collect_soiap.get_stress_step_soiap(
            stress_step_data, current_id, work_path+'log.strs')
        return stress_step_data
    elif rin.calc_code == 'LAMMPS':
        raise NotImplementedError('now only VASP')
    else:
        raise NotImplementedError('now only VASP')
