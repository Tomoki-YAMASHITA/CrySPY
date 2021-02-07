'''
Collect results in soiap
'''

import numpy as np

from . import structure as soiap_structure
from ... import utility
from ...IO import pkl_data
from ...IO import read_input as rin


def collect_soiap(current_id, work_path):
    # ---------- check optimization in current stage
    try:
        with open(work_path+rin.soiap_outfile, 'r') as fout:
            lines = fout.readlines()
        check_opt = 'not_yet'
        for i, line in enumerate(lines):
            if '*** QMD%loopc' in line:
                if ('QMD%frc converged.' in lines[i-2]
                        and 'QMD%strs converged.' in lines[i-1]):
                    check_opt = 'done'
                break
    except:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    magmom = np.nan    # magnetic moment is not calculated
    try:
        with open(work_path+'log.tote') as f:
            lines = f.readlines()
        energy = float(lines[-1].split()[2])    # in Hartree
        energy = energy * utility.hrt2ev        # Hartree to eV
        energy = energy/rin.natot               # eV/cell --> eV/atom
    except:
        energy = np.nan    # error
        print('    Structure ID {0}, could not obtain energy from {1}'.format(
            current_id, rin.soiap_outfile))

    # ---------- collect the last structure
    try:
        with open(work_path+'log.struc', 'r') as f:
            lines = f.readlines()
            lines = lines[-(rin.natot+5):]
        opt_struc = soiap_structure.from_file(lines)
    except:
        opt_struc = None

    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def get_energy_step_soiap(energy_step_data, current_id, filename):
    '''
    get energy step data in eV/atom

    energy_step_data[ID][stage][step]
    energy_step_data[ID][0] <-- stage 1
    energy_step_data[ID][1] <-- stage 2
    '''
    try:
        log_np = np.loadtxt(filename)
        energy_step = log_np[:, 4]    # collumn 4: Hartree/atom
        energy_step = energy_step * utility.hrt2ev    # eV/atom
    except:
        energy_step = None
        print('\n#### ID: {0}: failed to parse log.tote\n\n'.format(
            current_id))

    # ---------- append energy_step
    if energy_step_data.get(current_id) is None:
        energy_step_data[current_id] = []    # initialize
    energy_step_data[current_id].append(energy_step)

    # ---------- save energy_step_data
    pkl_data.save_energy_step(energy_step_data)

    # ---------- return
    return energy_step_data


def get_struc_step_soiap(struc_step_data, current_id, filename):
    '''
    get structure step data

    # ---------- args
    struc_step_data: (dict) the key is structure ID

    struc_step_data[ID][stage][step]
    struc_step_data[ID][0] <-- stage 1
    struc_step_data[ID][1] <-- stage 2
    '''
    # ---------- get struc step from log.struc
    try:
        # ------ read file
        with open(filename, 'r') as f:
            lines = f.readlines()
        # ------ init.
        struc_step = []
        # ------ loop for relaxation step
        tmp_lines = []
        for line in lines:
            tmp_lines.append(line)
            if len(tmp_lines) == rin.natot + 5:
                struc = soiap_structure.from_file(tmp_lines)
                struc_step.append(struc)
                tmp_lines = []    # clear
    except:
        struc_step = None
        print('\n#### ID: {0}: failed to parse log.struc\n\n'.format(
            current_id))

    # ---------- append struc_step
    if struc_step_data.get(current_id) is None:
        struc_step_data[current_id] = []    # initialize
    struc_step_data[current_id].append(struc_step)

    # ---------- save struc_step_data
    pkl_data.save_struc_step(struc_step_data)

    # ---------- return
    return struc_step_data


def get_force_step_soiap(force_step_data, current_id, filename):
    '''
    get force step data

    # ---------- args
    force_step_data: (dict) the key is structure ID

    force_step_data[ID][stage][step]
    force_step_data[ID][0] <-- stage 1
    force_step_data[ID][1] <-- stage 2
    '''
    # ---------- get force step from log.frc
    try:
        # ------ read file
        with open(filename, 'r') as f:
            lines = f.readlines()
        # ------ init
        force_step = []
        tmp_lines = []
        # ------ parse
        for line in lines:
            if 'forces' not in line:
                tmp_lines.append([float(x) for x in line.split()])
                if len(tmp_lines) == rin.natot:
                    tmp_lines = np.array(tmp_lines)
                    # Hartree/Bohr --> eV/ang
                    tmp_lines = tmp_lines * utility.hrt2ev / utility.bohr2ang
                    force_step.append(tmp_lines)
                    tmp_lines = []    # clear
    except:
        force_step = None
        print('\n#### ID: {0}: failed to parse log.frc\n\n'.format(
            current_id))

    # ---------- append force_step
    if force_step_data.get(current_id) is None:
        force_step_data[current_id] = []    # initialize
    force_step_data[current_id].append(force_step)

    # ---------- save force_step_data
    pkl_data.save_force_step(force_step_data)

    # ---------- return
    return force_step_data


def get_stress_step_soiap(stress_step_data, current_id, filename):
    '''
    get stress step data

    # ---------- args
    stress_step_data: (dict) the key is structure ID

    stress_step_data[ID][stage][step]
    stress_step_data[ID][0] <-- stage 1
    stress_step_data[ID][1] <-- stage 2
    '''
    # ---------- get stress step from log.strs
    try:
        # ------ read file
        with open(filename, 'r') as f:
            lines = f.readlines()
        # ------ init
        stress_step = []
        tmp_lines = []
        # ------ parse
        for line in lines:
            if 'QMD' not in line:
                tmp_lines.append([float(x) for x in line.split()])
                if len(tmp_lines) == 3:
                    tmp_lines = np.array(tmp_lines)
                    # Hartree/Bohr**3 --> eV/ang**3
                    tmp_lines = tmp_lines * utility.hrt2ev / utility.bohr2ang**3
                    stress_step.append(tmp_lines)
                    tmp_lines = []    # clear
    except:
        stress_step = None
        print('\n#### ID: {0}: failed to parse log.strs\n\n'.format(
            current_id))

    # ---------- append stress_step
    if stress_step_data.get(current_id) is None:
        stress_step_data[current_id] = []    # initialize
    stress_step_data[current_id].append(stress_step)

    # ---------- save stress_step_data
    pkl_data.save_stress_step(stress_step_data)

    # ---------- return
    return stress_step_data
