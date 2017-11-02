#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil

from pymatgen.io.vasp.sets import MITRelaxSet

from . import structure as qe_structure
from ...IO import out_kpts
from ...IO import pkl_data
from ...IO import read_input as rin


def next_stage_qe(stage, work_path, kpt_data, current_id):
    # ---------- skip_flag
    skip_flag = False

    # ---------- prepare QE files
    qe_files = [rin.qe_infile, rin.qe_outfile]
    for f in qe_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)

    # ---------- next structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(work_path+'prev_'+rin.qe_outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(work_path+'prev_'+rin.qe_infile)
        lines_atom = qe_structure.extract_atomic_positions(work_path+'prev_'+rin.qe_outfile)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(work_path+'prev_'+rin.qe_infile)
        structure = qe_structure.from_lines(lines_cell, lines_atom)
    except ValueError:
        skip_flag = True
        kpt_data[current_id].append(['skip'])
        pkl_data.save_kpt(kpt_data)
        out_kpts.write_kpts(kpt_data)
        print('    error in QE,  skip this structure')
        return skip_flag, kpt_data

    # ---------- copy the input file from ./calc_in
    finput = './calc_in/'+rin.qe_infile+'_{}'.format(stage)
    shutil.copyfile(finput, work_path+rin.qe_infile)

    # ---------- append structure info.
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kppvol[stage-1])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

    # ---------- kpt_data
    kpt_data[current_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts.write_kpts(kpt_data)

    # ---------- return
    return skip_flag, kpt_data


def next_struc_qe(init_struc_data, next_id, work_path, kpt_data):
    # ---------- copy files
    calc_inputs = [rin.qe_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.qe_infile else f
        if not os.path.isfile('./calc_in/' + ff):
            raise IOError('Could not find ./calc_in/' + ff)
        # ------ e.g. cp ./calc_in/xxxxx_1 work0001/xxxxx
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- append structure info. to the input file
    structure = init_struc_data[next_id]
    with open(work_path+rin.qe_infile, 'a') as fin:
        fin.write('\n')
    qe_structure.write(structure, work_path+rin.qe_infile, mode='a')

    # ---------- K_POINTS
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, rin.kppvol[0])
    with open(work_path+rin.qe_infile, 'a') as f:
        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' '.join(str(x) for x in kpoints.kpts[0]) + '  0 0 0\n')

    # ---------- kpt_data
    kpt_data[next_id] = []    # initialize
    kpt_data[next_id].append(kpoints.kpts[0])
    pkl_data.save_kpt(kpt_data)
    out_kpts.write_kpts(kpt_data)

    # ---------- return
    return kpt_data
