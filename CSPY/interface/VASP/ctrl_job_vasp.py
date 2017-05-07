#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from pymatgen import Structure
from pymatgen.io.vasp.sets import MITRelaxSet

from ...IO import read_input as rin


def next_stage_vasp(stage, work_path):
    #---------- prepare vasp files
    vasp_files = ['POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR']
    for f in vasp_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found ' +work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)
    shutil.copyfile(work_path+'prev_CONTCAR', work_path+'POSCAR')

    #---------- KPOINTS using pymatgen
    try:
        structure = Structure.from_file(work_path+'POSCAR')
    except ValueError:
        raise ValueError('Error in VASP. Check '+work_path+'. If you skip this structure, write "skip" in stat_job line 3.')
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kmesh[stage-1])
    kpoints.write_file(work_path+'KPOINTS')

    #---------- cp INCAR_? from ./calc_in
    fincar = './calc_in/INCAR_{}'.format(stage)
    shutil.copyfile(fincar, work_path+'INCAR')


def next_struc_vasp(init_struc_data, next_id, work_path):
    #---------- copy files
    calc_inputs = ['POTCAR', 'INCAR']
    for f in calc_inputs:
        if f == 'INCAR':
            ff = f+'_1'
        else:
            ff = f
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
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kmesh[0])
    kpoints.write_file(work_path+'KPOINTS')
