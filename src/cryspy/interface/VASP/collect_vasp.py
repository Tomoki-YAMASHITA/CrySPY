'''
Collect results in VASP
'''

from logging import getLogger
import os

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun, Outcar

from ...util import constants
from ...IO import pkl_data


logger = getLogger('cryspy')


def collect_vasp(cid, work_path, nat):
    # ---------- init
    opt_struc = None
    energy = np.nan
    magmom = np.nan
    check_opt = 'no_file'

    # ---------- read vasprun.xml
    try:
        vasprun = Vasprun(
            work_path + 'vasprun.xml',
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
    except Exception as e:
        logger.error(f'    Structure ID {cid}, failed to parse vasprun.xml: {e}')
        return opt_struc, energy, magmom, check_opt

    # ---------- check optimization
    if vasprun.converged:
        check_opt = 'done'
    else:
        check_opt = 'not_yet'

    # ---------- energy
    energy = vasprun.final_energy    # if failed, float('inf')
    if energy == float('inf'):
        logger.warning(f'    Structure ID {cid},'
              ' could not obtain energy from vasprun.xml')
        energy = np.nan
        return opt_struc, energy, magmom, check_opt
    energy = energy / sum(nat)    # eV/atom

    # ---------- magmom
    outcar = Outcar(work_path + 'OUTCAR')
    magmom = outcar.total_mag
    if magmom is None:
        magmom = np.nan

    # ---------- collect CONTCAR
    try:
        opt_struc = Structure.from_file(work_path + 'CONTCAR')
    except Exception:
        logger.error(f'    Structure ID {cid}, failed to parse CONTCAR')
        opt_struc = None
        energy = np.nan
        magmom = np.nan
        check_opt = 'no_file'

    # ---------- remove STOPCAR
    if os.path.isfile(work_path + 'STOPCAR'):
        os.remove(work_path + 'STOPCAR')

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def get_energy_step_vasp(energy_step_data, cid, work_path, nat):
    '''
    get energy step data in eV/atom

    energy_step_data[ID][stage][step]
    energy_step_data[ID][0] <-- stage 1
    energy_step_data[ID][1] <-- stage 2
    '''
    # ---------- read vasprun.xml
    try:
        vasprun = Vasprun(
            work_path + 'vasprun.xml',
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        energy_step = [step['e_0_energy'] for step in vasprun.ionic_steps]
        energy_step = np.array(energy_step)/sum(nat)
    except Exception as e:
        energy_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse in energy_step')

    # ---------- append energy_step
    if energy_step_data.get(cid) is None:
        energy_step_data[cid] = []    # initialize
    energy_step_data[cid].append(energy_step)

    # ---------- save energy_step_data
    pkl_data.save_energy_step(energy_step_data)

    # ---------- return
    return energy_step_data


def get_struc_step_vasp(struc_step_data, cid, work_path):
    '''
    get structure step data

    # ---------- args
    struc_step_data: (dict) the key is structure ID

    struc_step_data[ID][stage][step]
    struc_step_data[ID][0] <-- stage 1
    struc_step_data[ID][1] <-- stage 2
    '''
    # ---------- read vasprun.xml
    try:
        vasprun = Vasprun(
            work_path + 'vasprun.xml',
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        struc_step = vasprun.structures
    except Exception as e:
        struc_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse in struc_step')

    # ---------- append struc_step
    if struc_step_data.get(cid) is None:
        struc_step_data[cid] = []    # initialize
    struc_step_data[cid].append(struc_step)

    # ---------- save struc_step_data
    pkl_data.save_struc_step(struc_step_data)

    # ---------- return
    return struc_step_data


def get_force_step_vasp(force_step_data, cid, work_path):
    '''
    get force step data in eV/angstrom

    # ---------- args
    force_step_data: (dict) the key is structure ID

    force_step_data[ID][stage][step]
    force_step_data[ID][0] <-- stage 1
    force_step_data[ID][1] <-- stage 2
    '''
    # ---------- read vasprun.xml
    try:
        vasprun = Vasprun(
            work_path + 'vasprun.xml',
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        force_step = [step['forces'] for step in vasprun.ionic_steps]
    except Exception as e:
        force_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse in force_step')

    # ---------- append force_step
    if force_step_data.get(cid) is None:
        force_step_data[cid] = []    # initialize
    force_step_data[cid].append(force_step)

    # ---------- save force_step_data
    pkl_data.save_force_step(force_step_data)

    # ---------- return
    return force_step_data


def get_stress_step_vasp(stress_step_data, cid, work_path):
    '''
    get stress step data in eV/ang**3

    # ---------- args
    stress_step_data: (dict) the key is structure ID

    stress_step_data[ID][stage][step]
    stress_step_data[ID][0] <-- stage 1
    stress_step_data[ID][1] <-- stage 2
    '''
    # ---------- read vasprun.xml
    try:
        vasprun = Vasprun(
            work_path + 'vasprun.xml',
            parse_dos=False,
            parse_eigen=False,
            parse_potcar_file=False,
        )
        # kbar --> eV/ang**3
        stress_step = [step['stress'] * constants.KBAR2eV_ANG3 for step in vasprun.ionic_steps]
    except Exception as e:
        stress_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse in stress_step')

    # ---------- append stress_step
    if stress_step_data.get(cid) is None:
        stress_step_data[cid] = []    # initialize
    stress_step_data[cid].append(stress_step)

    # ---------- save stress_step_data
    pkl_data.save_stress_step(stress_step_data)

    # ---------- return
    return stress_step_data
