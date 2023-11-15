'''
Select an optimizer
    - VASP
    - QE
    - OMX
    - soiap
    - LAMMPS
    - ASE
'''

from logging import getLogger

from ..IO import read_input as rin

if rin.calc_code == 'VASP':
    from .VASP import calc_files_vasp, ctrl_job_vasp, collect_vasp
elif rin.calc_code == 'QE':
    from .QE import calc_files_qe, ctrl_job_qe, collect_qe
elif rin.calc_code == 'OMX':
    from .OMX import calc_files_OMX, ctrl_job_OMX, collect_OMX
elif rin.calc_code == 'soiap':
    from .soiap import calc_files_soiap, ctrl_job_soiap, collect_soiap
elif rin.calc_code == 'LAMMPS':
    from .LAMMPS import calc_files_lammps, ctrl_job_lammps, collect_lammps
elif rin.calc_code == 'ASE':
    from .ASE import calc_files_ase, ctrl_job_ase, collect_ase
elif rin.calc_code == 'ext':
    from .ext import collect_ext


logger = getLogger('cryspy')

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
    elif rin.calc_code == 'ASE':
        calc_files_ase.check_input_ase()
    elif rin.calc_code == 'ext':
        pass
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def next_stage(stage, work_path, nat, kpt_data=None, cid=None):
    if rin.calc_code == 'VASP':
        skip_flag, kpt_data = ctrl_job_vasp.next_stage_vasp(stage, work_path,
                                                            kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'QE':
        # qe needs the number of atoms
        skip_flag, kpt_data = ctrl_job_qe.next_stage_qe(stage, work_path,
                                                        nat, kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'OMX':
        # omx needs the number of atoms
        skip_flag, kpt_data = ctrl_job_OMX.next_stage_OMX(stage, work_path,
                                                          nat, kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'soiap':
        # soiap needs the number of atoms
        skip_flag = ctrl_job_soiap.next_stage_soiap(stage, work_path, nat)
        return skip_flag
    elif rin.calc_code == 'LAMMPS':
        # lammps needs the number of atoms
        skip_flag = ctrl_job_lammps.next_stage_lammps(stage, work_path, nat)
        return skip_flag
    elif rin.calc_code == 'ASE':
        skip_flag = ctrl_job_ase.next_stage_ase(stage, work_path)
        return skip_flag
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def collect(current_id, work_path, nat):
    if rin.calc_code == 'VASP':
        # vasp needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_vasp.collect_vasp(current_id, work_path, nat)
    elif rin.calc_code == 'QE':
        # qe needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_qe.collect_qe(current_id, work_path, nat)
    elif rin.calc_code == 'OMX':
        # OMX needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_OMX.collect_OMX(current_id, work_path, nat)
    elif rin.calc_code == 'soiap':
        # soiap needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_soiap.collect_soiap(current_id, work_path, nat)
    elif rin.calc_code == 'LAMMPS':
        # lammps needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_lammps.collect_lammps(current_id, work_path, nat)
    elif rin.calc_code == 'ASE':
        # ase needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_ase.collect_ase(current_id, work_path, nat)
    elif rin.calc_code == 'ext':
        ext_opt_struc_data, ext_energy_data = collect_ext.collect_ext()
        return ext_opt_struc_data, ext_energy_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def next_struc(structure, current_id, work_path, kpt_data=None):
    # args[0] <-- kpt_data
    if rin.calc_code == 'VASP':
        kpt_data = ctrl_job_vasp.next_struc_vasp(structure, current_id,
                                                 work_path, kpt_data)
        return kpt_data
    elif rin.calc_code == 'QE':
        kpt_data = ctrl_job_qe.next_struc_qe(structure, current_id,
                                             work_path, kpt_data)
        return kpt_data
    elif rin.calc_code == 'OMX':
        kpt_data = ctrl_job_OMX.next_struc_OMX(structure, current_id,
                                               work_path, kpt_data)
        return kpt_data
    elif rin.calc_code == 'soiap':
        ctrl_job_soiap.next_struc_soiap(structure, current_id, work_path)
    elif rin.calc_code == 'LAMMPS':
        ctrl_job_lammps.next_struc_lammps(structure, current_id, work_path)
    elif rin.calc_code == 'ASE':
        ctrl_job_ase.next_struc_ase(structure, current_id, work_path)
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_energy_step(energy_step_data, current_id, work_path, nat):
    if rin.calc_code == 'VASP':
        # vasp needs the number of atoms to collect results
        energy_step_data = collect_vasp.get_energy_step_vasp(
            energy_step_data, current_id, work_path, nat)
        return energy_step_data
    elif rin.calc_code == 'QE':
        # qe needs the number of atoms to collect results
        energy_step_data = collect_qe.get_energy_step_qe(energy_step_data,
                                                         current_id, work_path, nat)
        return energy_step_data
    elif rin.calc_code == 'soiap':
        energy_step_data = collect_soiap.get_energy_step_soiap(
            energy_step_data, current_id, work_path)
        return energy_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_struc_step(struc_step_data, current_id, work_path, nat):
    if rin.calc_code == 'VASP':
        struc_step_data = collect_vasp.get_struc_step_vasp(
            struc_step_data, current_id, work_path)
        return struc_step_data
    elif rin.calc_code == 'QE':
        # qe needs the number of atoms to collect results
        struc_step_data = collect_qe.get_struc_step_qe(
            struc_step_data, current_id, work_path, nat)
        return struc_step_data
    elif rin.calc_code == 'soiap':
        # soiap needs the number of atoms to collect results
        struc_step_data = collect_soiap.get_struc_step_soiap(
            struc_step_data, current_id, work_path, nat)
        return struc_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_force_step(force_step_data, current_id, work_path, nat):
    if rin.calc_code == 'VASP':
        force_step_data = collect_vasp.get_force_step_vasp(
            force_step_data, current_id, work_path)
        return force_step_data
    elif rin.calc_code == 'QE':
        # qe needs the number of atoms to collect results
        force_step_data = collect_qe.get_force_step_qe(
            force_step_data, current_id, work_path, nat)
        return force_step_data
    elif rin.calc_code == 'soiap':
        # soiap needs the number of atoms to collect results
        force_step_data = collect_soiap.get_force_step_soiap(
            force_step_data, current_id, work_path, nat)
        return force_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_stress_step(stress_step_data, current_id, work_path):
    if rin.calc_code == 'VASP':
        stress_step_data = collect_vasp.get_stress_step_vasp(
            stress_step_data, current_id, work_path)
        return stress_step_data
    elif rin.calc_code == 'QE':
        stress_step_data = collect_qe.get_stress_step_qe(
            stress_step_data, current_id, work_path)
        return stress_step_data
    elif rin.calc_code == 'soiap':
        stress_step_data = collect_soiap.get_stress_step_soiap(
            stress_step_data, current_id, work_path)
        return stress_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)
