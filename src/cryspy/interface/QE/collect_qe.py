from logging import getLogger

import numpy as np
from pymatgen.core import Structure

from . import structure as qe_structure
from ...util import constants
from ...IO import pkl_data


logger = getLogger('cryspy')

def collect_qe(rin, cid, work_path, nat):
    # ---------- check optimization in previous stage
    try:
        with open(work_path+rin.qe_outfile, 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'End final coordinates' in line:
                check_opt = 'done'
    except Exception as e:
        logger.warning(e)
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    natot = sum(nat)
    try:
        with open(work_path+rin.qe_outfile, 'r') as fpout:
            lines = fpout.readlines()
        energy = np.nan
        for line in reversed(lines):
            if rin.pv_term:
                condition = 'Final enthalpy' in line
            else:
                condition = line.startswith('!')
            if condition:
                energy = float(line.split()[-2])    # in Ry
                energy = energy * constants.RY2EV / float(natot)    # Ry/cell --> eV/atom
                break
        magmom = np.nan    # implemented by H. Sawahata 2020/10/04
        for line in reversed(lines):
            if line.find("total magnetization") >= 0:
                muB = line.split()
                magmom = float(muB[3])
                break
    except Exception as e:
        energy = np.nan    # error
        magmom = np.nan    # error
        logger.warning(f'{e}:    Structure ID {cid},'
                       f'could not obtain energy from {rin.qe_outfile}')

    # ---------- collect the last structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(
            work_path+rin.qe_outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(
                work_path+rin.qe_infile)
        lines_atom = qe_structure.extract_atomic_positions(
            work_path+rin.qe_outfile, nat)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(
                work_path+rin.qe_infile, nat)
        opt_struc = qe_structure.from_lines(lines_cell, lines_atom)

        # ------ opt_qe-structure
        with open('./data/opt_qe-structure', 'a') as fstruc:
            fstruc.write(f'# ID {cid}\n')
        qe_structure.write(opt_struc, './data/opt_qe-structure', mode='a')
    except Exception as e:
        logger.warning(e)
        opt_struc = None

    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan

    # ---------- return
    return opt_struc, energy, magmom, check_opt


def get_energy_step_qe(rin, energy_step_data, cid, work_path, nat):
    '''
    get energy step data in eV/atom

    energy_step_data[ID][stage][step]
    energy_step_data[ID][0] <-- stage 1
    energy_step_data[ID][1] <-- stage 2
    '''
    # ---------- natot
    natot = sum(nat)
    try:
        # ---------- read output file
        with open(work_path+rin.qe_outfile, 'r') as f:
            lines = f.readlines()
        # ---------- get energy step
        energy_step = []
        final_flag = False   # End final coordinates
        vc_flag = False      # vc-relax
        for line in lines:
            if line.startswith('!'):
                energy_step.append(line.split()[4])
            # ------ check opt and vc-relax
            if 'End final coordinates' in line:
                final_flag = True
            if 'CELL_PARAMETERS' in line:
                vc_flag = True
        # ------ delete last energy (after End final coordinates)
        if final_flag and vc_flag:
            energy_step.pop(-1)
        # ------ list --> array, Ry/cell --> eV/atom
        if not energy_step:
            energy_step = None    # if empty
            logger.warning(f'#### ID: {cid}: failed to parse energy_step')
        else:
            energy_step = constants.RY2EV / natot * np.array(energy_step,
                                                               dtype='float')
    except Exception as e:
        energy_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse energy_step')

    # ---------- append energy_step
    if energy_step_data.get(cid) is None:
        energy_step_data[cid] = []    # initialize
    energy_step_data[cid].append(energy_step)

    # ---------- save energy_step_data
    pkl_data.save_energy_step(energy_step_data)

    # ---------- return
    return energy_step_data


def get_struc_step_qe(rin, struc_step_data, cid, work_path, nat):
    '''
    get structure step data

    # ---------- args
    struc_step_data: (dict) the key is structure ID

    struc_step_data[ID][stage][step]
    struc_step_data[ID][0] <-- stage 1
    struc_step_data[ID][1] <-- stage 2
    '''
    try:
        struc_step = []
        # ------ init struc from pwscf.in
        _extract_struc_qe(work_path+rin.qe_infile, struc_step, nat)
        # ------ struc step from pwscf.out
        _extract_struc_qe(work_path+rin.qe_outfile, struc_step, nat)
        # ------ delete last structure due to duplication
        struc_step.pop(-1)
    except Exception as e:
        struc_step = None
        logger.warning(f'{e}:    ID {cid}: failed to parse in struc_step')

    # ---------- append struc_step_data
    if struc_step_data.get(cid) is None:
        struc_step_data[cid] = []    # initialize
    struc_step_data[cid].append(struc_step)

    # ---------- save struc_step_data
    pkl_data.save_struc_step(struc_step_data)

    # ---------- return
    return struc_step_data


def _extract_struc_qe(filename, struc_step, nat):
    # ---------- read a file
    with open(filename, 'r') as f:
        lines = f.readlines()
    # ---------- extract struc
    natot = sum(nat)
    read_cell = False
    read_coords = False
    vc_flag = False      # in case of vc-relax
    for line in lines:
        # ------ cell part
        if read_cell:
            lattice.append(line.split())
            if len(lattice) == 3:
                read_cell = False
                lattice = np.array(lattice, dtype='float')
        if 'CELL_PARAMETERS' in line:
            read_cell = True
            vc_flag = True
            lattice = []
        # ------ coords part
        if read_coords:
            lsplit = line.split()
            species.append(lsplit[0])
            coords.append(lsplit[1:])
            if len(coords) == natot:
                read_coords = False
                coords = np.array(coords, dtype='float')
                # ---- gen struc
                if not vc_flag:    # empty lattice, use init lattice
                    lattice = struc_step[0].lattice
                struc = Structure(lattice, species, coords)
                struc_step.append(struc)
        if 'ATOMIC_POSITIONS' in line:
            read_coords = True
            species = []
            coords = []


def get_force_step_qe(rin, force_step_data, cid, work_path, nat):
    '''
    get force step data in eV/angstrom

    # ---------- args
    force_step_data: (dict) the key is structure ID

    force_step_data[ID][stage][step]
    force_step_data[ID][0] <-- stage 1
    force_step_data[ID][1] <-- stage 2
    '''
    # ---------- natot
    natot = sum(nat)
    try:
        # ---------- read output file
        with open(work_path+rin.qe_outfile, 'r') as f:
            lines = f.readlines()
        # ---------- get force step
        force_step = []
        read_force = False
        final_flag = False   # End final coordinates
        vc_flag = False      # in case of vc-relax
        for line in lines:
            if 'atom    1 type' in line:
                read_force = True
                force = []
            if read_force:
                force.append(line.split()[6:])
                if len(force) == natot:
                    read_force = False
                    force_step.append(constants.RY2EV / constants.BOHR2ANG * np.array(
                        force, dtype='float'))
            # ------ check opt and vc-relax
            if 'End final coordinates' in line:
                final_flag = True
            if 'CELL_PARAMETERS' in line:
                vc_flag = True
        # ------ delete last energy (after End final coordinates)
        if final_flag and vc_flag:
            force_step.pop(-1)
        # ------ if empty
        if len(force_step) == 0:
            force_step = None
            logger.warning(f'#### ID: {cid}: failed to parse force_step')
    except Exception as e:
        force_step = None
        logger.warning(e + f'#### ID: {cid}: failed to parse in force_step')

    # ---------- append force_step
    if force_step_data.get(cid) is None:
        force_step_data[cid] = []    # initialize
    force_step_data[cid].append(force_step)

    # ---------- save force_step_data
    pkl_data.save_force_step(force_step_data)

    # ---------- return
    return force_step_data


def get_stress_step_qe(rin, stress_step_data, cid, work_path):
    '''
    get stress step data in eV/ang**3

    # ---------- args
    stress_step_data: (dict) the key is structure ID

    stress_step_data[ID][stage][step]
    stress_step_data[ID][0] <-- stage 1
    stress_step_data[ID][1] <-- stage 2
    '''
    try:
        # ---------- read output file
        with open(work_path+rin.qe_outfile, 'r') as f:
            lines = f.readlines()
        # ---------- get stress step
        stress_step = []
        read_stress = False
        final_flag = False   # End final coordinates
        vc_flag = False      # in case of vc-relax
        for line in lines:
            if read_stress:
                stress.append(line.split()[3:])
                if len(stress) == 3:
                    read_stress = False
                    stress_step.append(constants.KBAR2eV_ANG3 * np.array(
                        stress, dtype='float'))
            if 'total   stress  (Ry/bohr**3)' in line:
                read_stress = True
                stress = []
                # ------ check opt and vc-relax
            if 'End final coordinates' in line:
                final_flag = True
            if 'CELL_PARAMETERS' in line:
                vc_flag = True
        # ------ delete last energy (after End final coordinates)
        if final_flag and vc_flag:
            stress_step.pop(-1)
        # ------ if empty
        if len(stress_step) == 0:
            stress_step = None
            logger.warning(f'#### ID: {cid}: failed to parse stress_step')
    except Exception as e:
        stress_step = None
        logger.warning(e + f'#### ID: {cid}: failed to parse in stress_step')

    # ---------- append stress_step
    if stress_step_data.get(cid) is None:
        stress_step_data[cid] = []    # initialize
    stress_step_data[cid].append(stress_step)

    # ---------- save stress_step_data
    pkl_data.save_stress_step(stress_step_data)

    # ---------- return
    return stress_step_data
