'''
Collect results in soiap
'''

from logging import getLogger

import numpy as np

from . import structure as soiap_structure
from ...util import constants
from ...IO import pkl_data


logger = getLogger('cryspy')

def collect_soiap(rin, cid, work_path, nat):
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
    except Exception:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    magmom = np.nan    # magnetic moment is not calculated
    try:
        with open(work_path+'log.tote') as f:
            lines = f.readlines()
        energy = float(lines[-1].split()[4])    # in Hartree/atom
        energy = energy * constants.HRT2EV        # Hartree/atom to eV/atom
    except Exception as e:
        energy = np.nan    # error
        logger.warning(f'{e}:    Structure ID {cid},'
                    f' could not obtain energy from {rin.soiap_outfile}')

    # ---------- collect the last structure
    natot = sum(nat)
    try:
        with open(work_path+'log.struc', 'r') as f:
            lines = f.readlines()
            lines = lines[-(natot+5):]
        opt_struc = soiap_structure.from_file(rin, lines, nat)
    except Exception:
        opt_struc = None

    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def get_energy_step_soiap(energy_step_data, cid, work_path):
    '''
    get energy step data in eV/atom

    energy_step_data[ID][stage][step]
    energy_step_data[ID][0] <-- stage 1
    energy_step_data[ID][1] <-- stage 2

    In soiap, collect energy step data only when loopa == 1.
        This is because other data (struc, force, stress)
        are output only when loopa == 1
        see, https://github.com/nbsato/soiap/blob/master/doc/instructions.md
    '''
    try:
        energy_step = []
        with open(work_path+'log.tote') as f:
            lines = f.readlines()
        for line in lines:
            if line.split()[1] == '1':
                energy_step.append(line.split()[4])    # collumn 4: Hartree/atom
        energy_step = np.array(energy_step, dtype='float') * constants.HRT2EV
    except Exception as e:
        energy_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse log.tote')

    # ---------- append energy_step
    if energy_step_data.get(cid) is None:
        energy_step_data[cid] = []    # initialize
    energy_step_data[cid].append(energy_step)

    # ---------- save energy_step_data
    pkl_data.save_energy_step(energy_step_data)

    # ---------- return
    return energy_step_data


def get_struc_step_soiap(rin, struc_step_data, cid, work_path, nat):
    '''
    get structure step data

    # ---------- args
    struc_step_data: (dict) the key is structure ID

    struc_step_data[ID][stage][step]
    struc_step_data[ID][0] <-- stage 1
    struc_step_data[ID][1] <-- stage 2
    '''
    # ---------- get struc step from log.struc
    natot = sum(nat)
    try:
        # ------ read file
        with open(work_path+'log.struc', 'r') as f:
            lines = f.readlines()
        # ------ init.
        struc_step = []
        # ------ loop for relaxation step
        tmp_lines = []
        for line in lines:
            tmp_lines.append(line)
            if len(tmp_lines) == natot + 5:
                struc = soiap_structure.from_file(rin, tmp_lines, nat)
                struc_step.append(struc)
                tmp_lines = []    # clear
    except Exception as e:
        struc_step = None
        logger.warning(f'{e}    ID {cid}: failed to parse log.struc')

    # ---------- append struc_step
    if struc_step_data.get(cid) is None:
        struc_step_data[cid] = []    # initialize
    struc_step_data[cid].append(struc_step)

    # ---------- save struc_step_data
    pkl_data.save_struc_step(struc_step_data)

    # ---------- return
    return struc_step_data


def get_force_step_soiap(force_step_data, cid, work_path, nat):
    '''
    get force step data in eV/angstrom

    # ---------- args
    force_step_data: (dict) the key is structure ID

    force_step_data[ID][stage][step]
    force_step_data[ID][0] <-- stage 1
    force_step_data[ID][1] <-- stage 2
    '''
    # ---------- get force step from log.frc
    natot = sum(nat)
    try:
        # ------ read file
        with open(work_path+'log.frc', 'r') as f:
            lines = f.readlines()
        # ------ init
        force_step = []
        tmp_lines = []
        # ------ parse
        for line in lines:
            if 'forces' not in line:
                tmp_lines.append([float(x) for x in line.split()])
                if len(tmp_lines) == natot:
                    tmp_lines = np.array(tmp_lines)
                    # Hartree/Bohr --> eV/ang
                    tmp_lines = tmp_lines * constants.HRT2EV / constants.BOHR2ANG
                    force_step.append(tmp_lines)
                    tmp_lines = []    # clear
    except Exception as e:
        force_step = None
        logger.warning(f'{e}    ID {cid}: failed to parse log.frc')

    # ---------- append force_step
    if force_step_data.get(cid) is None:
        force_step_data[cid] = []    # initialize
    force_step_data[cid].append(force_step)

    # ---------- save force_step_data
    pkl_data.save_force_step(force_step_data)

    # ---------- return
    return force_step_data


def get_stress_step_soiap(stress_step_data, cid, work_path):
    '''
    get stress step data in eV/ang**3

    # ---------- args
    stress_step_data: (dict) the key is structure ID

    stress_step_data[ID][stage][step]
    stress_step_data[ID][0] <-- stage 1
    stress_step_data[ID][1] <-- stage 2
    '''
    # ---------- get stress step from log.strs
    try:
        # ------ read file
        with open(work_path+'log.strs', 'r') as f:
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
                    tmp_lines = tmp_lines * constants.HRT2EV / constants.BOHR2ANG**3
                    stress_step.append(tmp_lines)
                    tmp_lines = []    # clear
    except Exception as e:
        stress_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse log.strs')

    # ---------- append stress_step
    if stress_step_data.get(cid) is None:
        stress_step_data[cid] = []    # initialize
    stress_step_data[cid].append(stress_step)

    # ---------- save stress_step_data
    pkl_data.save_stress_step(stress_step_data)

    # ---------- return
    return stress_step_data
