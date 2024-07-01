from logging import getLogger

# ---------- import later
# from .VASP import calc_files_vasp, ctrl_job_vasp, collect_vasp
# from .QE import calc_files_qe, ctrl_job_qe, collect_qe
# from .OMX import calc_files_OMX, ctrl_job_OMX, collect_OMX
# from .soiap import calc_files_soiap, ctrl_job_soiap, collect_soiap
# from .LAMMPS import calc_files_lammps, ctrl_job_lammps, collect_lammps
# from .ASE import calc_files_ase, ctrl_job_ase, collect_ase


logger = getLogger('cryspy')


def check_calc_files(rin):
    if rin.calc_code == 'VASP':
        from .VASP import calc_files_vasp
        calc_files_vasp.check_input_vasp(rin)
    elif rin.calc_code == 'QE':
        from .QE import calc_files_qe
        calc_files_qe.check_input_qe(rin)
    elif rin.calc_code == 'OMX':
        from .OMX import calc_files_OMX
        calc_files_OMX.check_input_OMX(rin)
    elif rin.calc_code == 'soiap':
        from .soiap import calc_files_soiap
        calc_files_soiap.check_input_soiap(rin)
    elif rin.calc_code == 'LAMMPS':
        from .LAMMPS import calc_files_lammps
        calc_files_lammps.check_input_lammps(rin)
    elif rin.calc_code == 'ASE':
        from .ASE import calc_files_ase
        calc_files_ase.check_input_ase(rin)
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def next_stage(rin, stage, work_path, nat, kpt_data=None, cid=None):
    if rin.calc_code == 'VASP':
        from .VASP import ctrl_job_vasp
        skip_flag, kpt_data = ctrl_job_vasp.next_stage_vasp(rin, stage, work_path,
                                                            kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'QE':
        from .QE import ctrl_job_qe
        # qe needs the number of atoms
        skip_flag, kpt_data = ctrl_job_qe.next_stage_qe(rin, stage, work_path,
                                                        nat, kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'OMX':
        from .OMX import ctrl_job_OMX
        # omx needs the number of atoms
        skip_flag, kpt_data = ctrl_job_OMX.next_stage_OMX(rin, stage, work_path,
                                                          nat, kpt_data, cid)
        return skip_flag, kpt_data
    elif rin.calc_code == 'soiap':
        from .soiap import ctrl_job_soiap
        # soiap needs the number of atoms
        skip_flag = ctrl_job_soiap.next_stage_soiap(rin, stage, work_path, nat)
        return skip_flag
    elif rin.calc_code == 'LAMMPS':
        from .LAMMPS import ctrl_job_lammps
        # lammps needs the number of atoms
        skip_flag = ctrl_job_lammps.next_stage_lammps(rin, stage, work_path, nat)
        return skip_flag
    elif rin.calc_code == 'ASE':
        from .ASE import ctrl_job_ase
        skip_flag = ctrl_job_ase.next_stage_ase(rin, stage, work_path)
        return skip_flag
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def next_struc(rin, structure, cid, work_path, nat, kpt_data=None):
    # args[0] <-- kpt_data
    if rin.calc_code == 'VASP':
        from .VASP import ctrl_job_vasp
        kpt_data = ctrl_job_vasp.next_struc_vasp(rin, structure, cid,
                                                 work_path, kpt_data)
        return kpt_data
    elif rin.calc_code == 'QE':
        from .QE import ctrl_job_qe
        # qe needs the number of atoms
        kpt_data = ctrl_job_qe.next_struc_qe(rin, structure, cid,
                                             work_path, nat, kpt_data)
        return kpt_data
    elif rin.calc_code == 'OMX':
        from .OMX import ctrl_job_OMX
        kpt_data = ctrl_job_OMX.next_struc_OMX(rin, structure, cid,
                                               work_path, nat, kpt_data)
        return kpt_data
    elif rin.calc_code == 'soiap':
        from .soiap import ctrl_job_soiap
        ctrl_job_soiap.next_struc_soiap(rin, structure, cid, work_path)
    elif rin.calc_code == 'LAMMPS':
        from .LAMMPS import ctrl_job_lammps
        # lammps needs the number of atoms
        ctrl_job_lammps.next_struc_lammps(rin, structure, cid, work_path ,nat)
    elif rin.calc_code == 'ASE':
        from .ASE import ctrl_job_ase
        ctrl_job_ase.next_struc_ase(rin, structure, cid, work_path)
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def collect(rin, cid, work_path, nat):
    if rin.calc_code == 'VASP':
        from .VASP import collect_vasp
        # vasp needs the number of atoms to collect results
        # rin is not used in collect_vasp
        opt_struc, energy, magmom, check_opt = \
            collect_vasp.collect_vasp(cid, work_path, nat)
    elif rin.calc_code == 'QE':
        from .QE import collect_qe
        # qe needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_qe.collect_qe(rin, cid, work_path, nat)
    elif rin.calc_code == 'OMX':
        from .OMX import collect_OMX
        # OMX needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_OMX.collect_OMX(rin, cid, work_path, nat)
    elif rin.calc_code == 'soiap':
        from .soiap import collect_soiap
        # soiap needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_soiap.collect_soiap(rin, cid, work_path, nat)
    elif rin.calc_code == 'LAMMPS':
        from .LAMMPS import collect_lammps
        # lammps needs the number of atoms to collect results
        opt_struc, energy, magmom, check_opt = \
            collect_lammps.collect_lammps(rin, cid, work_path, nat)
    elif rin.calc_code == 'ASE':
        from .ASE import collect_ase
        # ase needs the number of atoms to collect results
        # rin is not used in collect_ase
        opt_struc, energy, magmom, check_opt = \
            collect_ase.collect_ase(cid, work_path, nat)
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def get_energy_step(rin, energy_step_data, cid, work_path, nat):
    from .VASP import collect_vasp
    if rin.calc_code == 'VASP':
        # vasp needs the number of atoms to collect results
        # rin is not used in collect_vasp
        energy_step_data = collect_vasp.get_energy_step_vasp(
            energy_step_data, cid, work_path, nat)
        return energy_step_data
    elif rin.calc_code == 'QE':
        from .QE import collect_qe
        # qe needs the number of atoms to collect results
        energy_step_data = collect_qe.get_energy_step_qe(rin, energy_step_data,
                                                         cid, work_path, nat)
        return energy_step_data
    elif rin.calc_code == 'soiap':
        from .soiap import collect_soiap
        # rin is not used in collect_soiap.get_energy_step_soiap
        energy_step_data = collect_soiap.get_energy_step_soiap(
            energy_step_data, cid, work_path)
        return energy_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_struc_step(rin, struc_step_data, cid, work_path, nat):
    if rin.calc_code == 'VASP':
        from .VASP import collect_vasp
        # rin is not used in collect_vasp
        struc_step_data = collect_vasp.get_struc_step_vasp(
            struc_step_data, cid, work_path)
        return struc_step_data
    elif rin.calc_code == 'QE':
        from .QE import collect_qe
        # qe needs the number of atoms to collect results
        struc_step_data = collect_qe.get_struc_step_qe(
            rin, struc_step_data, cid, work_path, nat)
        return struc_step_data
    elif rin.calc_code == 'soiap':
        from .soiap import collect_soiap
        # soiap needs the number of atoms to collect results
        struc_step_data = collect_soiap.get_struc_step_soiap(
            rin, struc_step_data, cid, work_path, nat)
        return struc_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_force_step(rin, force_step_data, cid, work_path, nat):
    if rin.calc_code == 'VASP':
        from .VASP import collect_vasp
        # rin is not used in collect_vasp
        force_step_data = collect_vasp.get_force_step_vasp(
            force_step_data, cid, work_path)
        return force_step_data
    elif rin.calc_code == 'QE':
        from .QE import collect_qe
        # qe needs the number of atoms to collect results
        force_step_data = collect_qe.get_force_step_qe(
            rin, force_step_data, cid, work_path, nat)
        return force_step_data
    elif rin.calc_code == 'soiap':
        from .soiap import collect_soiap
        # soiap needs the number of atoms to collect results
        # rin is not used in collect_soiap.get_force_step_soiap
        force_step_data = collect_soiap.get_force_step_soiap(
            force_step_data, cid, work_path, nat)
        return force_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)


def get_stress_step(rin, stress_step_data, cid, work_path):
    if rin.calc_code == 'VASP':
        from .VASP import collect_vasp
        # rin is not used in collect_vasp
        stress_step_data = collect_vasp.get_stress_step_vasp(
            stress_step_data, cid, work_path)
        return stress_step_data
    elif rin.calc_code == 'QE':
        from .QE import collect_qe
        stress_step_data = collect_qe.get_stress_step_qe(
            rin, stress_step_data, cid, work_path)
        return stress_step_data
    elif rin.calc_code == 'soiap':
        from .soiap import collect_soiap
        # rin is not used in collect_soiap.get_stress_step_soiap
        stress_step_data = collect_soiap.get_stress_step_soiap(
            stress_step_data, cid, work_path)
        return stress_step_data
    else:
        logger.error(f'{rin.calc_code}: not implemented yet')
        raise SystemExit(1)
