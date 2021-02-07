'''
Collect results in VASP
'''

import os
import xml.etree.ElementTree as ET

import numpy as np
from pymatgen import Structure

from ... import utility
from ...IO import pkl_data
from ...IO import read_input as rin


def collect_vasp(current_id, work_path):
    # ---------- check optimization
    check_opt = check_opt_vasp(work_path+'OUTCAR')
    # ---------- obtain energy and magmom
    energy, magmom = get_energy_magmom_vasp(work_path)
    if np.isnan(energy):
        print('    Structure ID {0},'
              ' could not obtain energy from OSZICAR'.format(current_id))
    # ---------- collect CONTCAR
    opt_struc = get_opt_struc_vasp(work_path+'CONTCAR')
    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan
    # ---------- remove STOPCAR
    if os.path.isfile(work_path+'STOPCAR'):
        os.remove(work_path+'STOPCAR')
    # ---------- return
    return opt_struc, energy, magmom, check_opt


def check_opt_vasp(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'reached required accuracy' in line:
                check_opt = 'done'
    except:
        check_opt = 'no_file'
    return check_opt


def get_energy_magmom_vasp(work_path):
    # ---------- obtain energy and magmom
    energy = np.nan
    magmom = np.nan
    try:
        with open(work_path+'OSZICAR', 'r') as foszi:
            oszi = foszi.readlines()
        if 'F=' in oszi[-1]:
            energy = float(oszi[-1].split()[2])    # free energy (eV/cell)
            energy = energy/float(rin.natot)       # eV/atom
            if 'mag=' in oszi[-1]:
                magmom = float(oszi[-1].split()[-1])    # total magnetic moment
    except:
        pass
    # ---------- return
    return energy, magmom


def get_opt_struc_vasp(file_name):
    try:
        opt_struc = Structure.from_file(file_name)
    except:
        opt_struc = None
    return opt_struc


def get_energy_step_vasp(energy_step_data, current_id, filename):
    '''
    get energy step data in eV/atom
    '''
    # ---------- get energy step from vasprun
    try:
        # ------ read file
        tree = ET.parse(filename)
        root = tree.getroot()
        # ------ children nodes: calculation
        cals = root.findall('calculation')
        # ------ init.
        energy_step = []
        # ------ loop for relaxation step
        for cal in cals:
            eng = cal.find('energy')    # first 'energy' child node
            fr_eng = eng.find('i')    # first 'i' tag is free energy

            if fr_eng.attrib['name'] == 'e_fr_energy':
                energy_step.append(fr_eng.text)
            else:
                raise ValueError('bug')
        # ------ list, str --> array
        energy_step = np.array(energy_step, dtype='float')/float(rin.natot)
    except:
        energy_step = None
        print('\n#### ID: {0}: failed to parse vasprun.xml\n\n'.format(
            current_id))

    # ---------- append energy_step
    if energy_step_data.get(current_id) is None:
        energy_step_data[current_id] = []    # initialize
    energy_step_data[current_id].append(energy_step)

    # ---------- save energy_step_data
    pkl_data.save_energy_step(energy_step_data)

    # ---------- return
    return energy_step_data


def get_struc_step_vasp(struc_step_data, current_id, filename):
    # ---------- get struc step from vasprun
    try:
        # ------ read file
        tree = ET.parse(filename)
        root = tree.getroot()
        # ------ get atom list
        atoms = root.findall("atominfo/array[@name='atoms']/set/rc")
        atomlist = []
        for atom in atoms:
            atomlist.append(atom.find('c').text)
        # ------ children nodes: calculation
        cals = root.findall('calculation')
        # ------ init.
        struc_step = []
        # ------ loop for relaxation step
        for cal in cals:
            # -- lattice
            basis = cal.findall("structure/crystal/varray[@name='basis']/v")
            lattice = []
            for a in basis:
                lattice.append([float(x) for x in a.text.split()])
            # -- positions
            positions = cal.findall("structure/varray/[@name='positions']/v")
            incoord = []
            for a in positions:
                incoord.append([float(x) for x in a.text.split()])
            # -- structure in pymatgen format
            struc = Structure(lattice, atomlist, incoord)
            # -- append
            struc_step.append(struc)
    except:
        struc_step = None
        print('\n#### ID: {0}: failed to parse vasprun.xml\n\n'.format(
            current_id))

    # ---------- append struc_step
    if struc_step_data.get(current_id) is None:
        struc_step_data[current_id] = []    # initialize
    struc_step_data[current_id].append(struc_step)

    # ---------- save struc_step_data
    pkl_data.save_struc_step(struc_step_data)

    # ---------- return
    return struc_step_data


def get_force_step_vasp(force_step_data, current_id, filename):
    # ---------- get force step from vasprun
    try:
        # ------ read file
        tree = ET.parse(filename)
        root = tree.getroot()
        # ------ children nodes: calculation
        cals = root.findall('calculation')
        # ------ init.
        force_step = []
        # ------ loop for ralaxation step
        for cal in cals:
            varrays = cal.findall('varray')
            # -- init
            force = []
            # -- varrays[0]: force, varrays[1]: stress
            for varray in varrays:
                vs = varray.findall('v')
                # loop for v
                for v in vs:
                    if varray.attrib['name'] == 'forces':
                        force.append(v.text.split())
            # -- list, str --> array
            force = np.array(force, dtype='float')
            # -- appned force_step
            force_step.append(force)
    except:
        force_step = None
        print('\n#### ID: {0}: failed to parse vasprun.xml\n\n'.format(
            current_id))

    # ---------- append force_step
    if force_step_data.get(current_id) is None:
        force_step_data[current_id] = []    # initialize
    force_step_data[current_id].append(force_step)

    # ---------- save force_step_data
    pkl_data.save_force_step(force_step_data)

    # ---------- return
    return force_step_data


def get_stress_step_vasp(stress_step_data, current_id, filename):
    # ---------- get stress step from vasprun
    try:
        # ------ read file
        tree = ET.parse(filename)
        root = tree.getroot()
        # ------ children nodes: calculation
        cals = root.findall('calculation')
        # ------ init.
        stress_step = []
        # ------ loop for ralaxation step
        for cal in cals:
            varrays = cal.findall('varray')
            # -- init
            stress = []
            # -- varrays[0]: force, varrays[1]: stress
            for varray in varrays:
                vs = varray.findall('v')
                # loop for v
                for v in vs:
                    if varray.attrib['name'] == 'stress':
                        stress.append(v.text.split())
            # -- list, str --> array
            stress = np.array(stress, dtype='float')
            # -- kbar --> eV/ang**3
            stress = stress * utility.kbar2ev_ang3
            # -- appned stress_step
            stress_step.append(stress)
    except:
        stress_step = None
        print('\n#### ID: {0}: failed to parse vasprun.xml\n\n'.format(
            current_id))

    # ---------- append stress_step
    if stress_step_data.get(current_id) is None:
        stress_step_data[current_id] = []    # initialize
    stress_step_data[current_id].append(stress_step)

    # ---------- save stress_step_data
    pkl_data.save_stress_step(stress_step_data)

    # ---------- return
    return stress_step_data
