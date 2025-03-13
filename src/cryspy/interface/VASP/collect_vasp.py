'''
Collect results in VASP
'''

from logging import getLogger
import os
import xml.etree.ElementTree as ET

import numpy as np
from pymatgen.core import Structure

from ...util import constants
from ...IO import pkl_data


logger = getLogger('cryspy')

def collect_vasp(cid, work_path, nat):
    # ---------- check optimization
    check_opt = check_opt_vasp(work_path+'OUTCAR')
    # ---------- obtain energy and magmom
    energy, magmom = get_energy_magmom_vasp(work_path, nat)
    if np.isnan(energy):
        logger.warning(f'    Structure ID {cid},'
              ' could not obtain energy from OSZICAR')
    # ---------- collect CONTCAR
    try:
        opt_struc = Structure.from_file(work_path+'CONTCAR')
    except Exception:
        opt_struc = None
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
    except Exception:
        check_opt = 'no_file'
    return check_opt


def get_energy_magmom_vasp(work_path, nat):
    # ---------- natot
    natot = sum(nat)
    # ---------- obtain energy and magmom
    energy = np.nan
    magmom = np.nan
    try:
        with open(work_path+'OSZICAR', 'r') as foszi:
            oszi = foszi.readlines()
        if 'F=' in oszi[-1]:
            energy = float(oszi[-1].split()[2])    # free energy (eV/cell)
            energy = energy/float(natot)       # eV/atom
            if 'mag=' in oszi[-1]:
                magmom = float(oszi[-1].split()[-1])    # total magnetic moment
    except Exception:
        pass
    # ---------- return
    return energy, magmom


def get_energy_step_vasp(energy_step_data, cid, work_path, nat):
    '''
    get energy step data in eV/atom

    energy_step_data[ID][stage][step]
    energy_step_data[ID][0] <-- stage 1
    energy_step_data[ID][1] <-- stage 2
    '''
    # ---------- natot
    natot = sum(nat)
    # ---------- get energy step from vasprun
    try:
        # ------ read file
        tree = ET.parse(work_path+'vasprun.xml')
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
                logger.error('bug in get_energy_step_vasp')
                raise SystemExit(1)
        # ------ list, str --> array
        energy_step = np.array(energy_step, dtype='float')/float(natot)
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
    # ---------- get struc step from vasprun
    try:
        # ------ read file
        tree = ET.parse(work_path+'vasprun.xml')
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
    # ---------- get force step from vasprun
    try:
        # ------ read file
        tree = ET.parse(work_path+'vasprun.xml')
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
    # ---------- get stress step from vasprun
    try:
        # ------ read file
        tree = ET.parse(work_path+'vasprun.xml')
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
            stress = stress * constants.KBAR2eV_ANG3
            # -- appned stress_step
            stress_step.append(stress)
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
