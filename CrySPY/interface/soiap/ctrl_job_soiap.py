#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import shutil

from . import structure as soiap_structure
from ...IO import read_input as rin


def next_stage_soiap(stage, work_path):
    # ---------- skip_flag
    skip_flag = False

    # ---------- prepare soiap files
    soiap_files = [rin.soiap_infile, rin.soiap_outfile, rin.soiap_cif,
                   'log.struc', 'log.tote', 'log.frc', 'log.strs']
    for f in soiap_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)

    # ---------- copy the input file from ./calc_in
    finfile = './calc_in/'+rin.soiap_infile+'_{}'.format(stage)
    shutil.copyfile(finfile, work_path+rin.soiap_infile)

    # ---------- generate the CIF file
    try:
        structure = soiap_structure.from_file(work_path+'prev_log.struc')
    except ValueError:
        skip_flag = True
        print('    error in soiap,  skip this structure')
        return skip_flag
    with open(work_path+'prev_'+rin.soiap_cif, 'r') as f:
        lines = f.readlines()
    title = lines[0][5:]    # string following 'data_'
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symtoleR,
                          title=title)

    # ---------- return
    return skip_flag


def next_struc_soiap(init_struc_data, next_id, work_path):
    # ---------- copy files
    calc_inputs = [rin.soiap_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.soiap_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            raise IOError('Could not find ./calc_in/'+ff)
        # ------ e.g. cp ./calc_in/xxxxx_1 work0001/xxxxx
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate the CIF file
    structure = init_struc_data[next_id]
    soiap_structure.write(structure,
                          work_path+rin.soiap_cif,
                          symprec=rin.symtoleI,
                          title='ID_{0:d}'.format(next_id))
