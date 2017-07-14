#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from pymatgen import Structure
from pymatgen.io.vasp.sets import MITRelaxSet

from ...IO import out_kpts
from ...IO import pkl_data
from ...IO import read_input as rin


def next_stage_vasp(stage, work_path, kpt_data, current_id):
    #---------- skip_flag
    skip_flag = False

    #---------- prepare vasp files
    vasp_files = ['POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR']
    for f in vasp_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)
    shutil.copyfile(work_path+'prev_CONTCAR', work_path+'POSCAR')

    #---------- remove STOPCAR
    if os.path.isfile(work_path+'STOPCAR'):
        os.remove(work_path+'STOPCAR')

    #---------- KPOINTS using pymatgen
    try:
        structure = Structure.from_file(work_path+'POSCAR')
    except ValueError:
        skip_flag = True
        kpt_data[current_id].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts.write_kpts(kpt_data)
        print('    error in VASP,  skip this structure')
        return skip_flag, kpt_data
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[stage-1],
                                                           rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    #---------- kpt_data
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts.write_kpts(kpt_data)

    #---------- cp INCAR_? from ./calc_in
    fincar = './calc_in/INCAR_{}'.format(stage)
    shutil.copyfile(fincar, work_path+'INCAR')

    #---------- return
    return skip_flag, kpt_data


def next_struc_vasp(init_struc_data, next_id, work_path, kpt_data):
    #---------- copy files
    calc_inputs = ['POTCAR', 'INCAR']
    for f in calc_inputs:
        ff = f+'_1' if f == 'INCAR' else f
        if not os.path.isfile('./calc_in/' + ff):
            raise IOError('Could not find ./calc_in/' + ff)
        #----- e.g. cp ./calc_in/INCAR_1 work0001/INCAR
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    #---------- generate POSCAR
    structure = init_struc_data[next_id]
    structure.to(fmt='poscar', filename=work_path+'POSCAR')
    if not os.path.isfile(work_path+'POSCAR'):
        raise IOError('Could not find {}POSCAR'.format(work_path))

    #---------- Change the title of POSCAR
    with open(work_path+'POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = 'ID_{}\n'.format(next_id)
    with open(work_path+'POSCAR', 'w') as f:
        for line in lines:
            f.write(line)

    #---------- generate KPOINTS using pymatgen
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure,
                                                           rin.kppvol[0],
                                                           rin.force_gamma)
    kpoints.write_file(work_path+'KPOINTS')

    #---------- kpt_data
    kpt_data[next_id] = []    # initialize
    kpt_data[next_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts.write_kpts(kpt_data)

    #---------- return
    return kpt_data
