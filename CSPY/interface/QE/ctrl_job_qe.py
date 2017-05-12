#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from pymatgen.io.vasp.sets import MITRelaxSet

from . import structure as qe_structure
from ...IO import read_input as rin


def next_stage_qe(stage, work_path):
    #---------- prepare QE files
    qe_files = [rin.infile, rin.outfile]
    for f in qe_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)

    #---------- next structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(work_path+'prev_'+rin.outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(work_path+'prev_'+rin.infile)
        lines_atom = qe_structure.extract_atomic_positions(work_path+'prev_'+rin.outfile)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(work_path+'prev_'+rin.infile)
        structure = qe_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        raise ValueError('Error in QE. Check '+work_path+'. If you skip this structure, write "skip" in stat_job line 3.')

    #---------- copy the input file from ./calc_in
    finput = './calc_in/'+rin.infile+'_{}'.format(stage)
    shutil.copyfile(finput, work_path+rin.infile)

    #---------- append structure info.
    with open(work_path+rin.infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.infile, mode='a')

    #---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kmesh[stage-1])
    with open(work_path+rin.infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')


def next_struc_qe(init_struc_data, next_id, work_path):
    #---------- copy files
    calc_inputs = [rin.infile] + rin.pseudopots
    for f in calc_inputs:
        if f == rin.infile:
            ff = f+'_1'
        else:
            ff = f
        if not os.path.isfile('./calc_in/' + ff):
            raise IOError('Could not find ./calc_in/' + ff)
        #----- e.g. cp ./calc_in/xxxxx_1 work0001/xxxxx
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    #---------- append structure info. to the input file
    structure = init_struc_data[next_id]
    with open(work_path+rin.infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.infile, mode='a')

    #---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kmesh[0])
    with open(work_path+rin.infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')
