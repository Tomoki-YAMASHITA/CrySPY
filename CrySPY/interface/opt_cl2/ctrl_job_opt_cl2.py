#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import shutil

from . import structure as opt_cl2_structure
from ...IO import read_input as rin


def next_stage_opt_cl2(stage, work_path):
    #---------- skip_flag
    skip_flag = False

    #---------- prepare opt_cl2 files
    opt_cl2_files = [rin.opt_cl2_infile, rin.opt_cl2_outfile, rin.opt_cl2_cif, 'log.struc', 'log.tote', 'log.frc', 'log.strs']
    for f in opt_cl2_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'prev_'+f)

    #---------- copy the input file from ./calc_in
    finfile = './calc_in/'+rin.opt_cl2_infile+'_{}'.format(stage)
    shutil.copyfile(finfile, work_path+rin.opt_cl2_infile)

    #---------- generate the CIF file
    try:
        structure = opt_cl2_structure.from_file(work_path+'prev_log.struc')
    except ValueError:
        skip_flag = True
        print('    error in opt_cl2,  skip this structure')
        return skip_flag
    with open(work_path+'prev_'+rin.opt_cl2_cif, 'r') as f:
        lines = f.readlines()
    title = lines[0][5:]    # string following 'data_'
    opt_cl2_structure.write(structure,
                            work_path+rin.opt_cl2_cif,
                            symprec=rin.symtoleR,
                            title=title)

    #---------- return
    return skip_flag


def next_struc_opt_cl2(init_struc_data, next_id, work_path):
    #---------- copy files
    calc_inputs = [rin.opt_cl2_infile]
    for f in calc_inputs:
        ff = f+'_1' if f == rin.opt_cl2_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            raise IOError('Could not find ./calc_in/'+ff)
        #----- e.g. cp ./calc_in/xxxxx_1 work0001/xxxxx
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    #---------- generate the CIF file
    structure = init_struc_data[next_id]
    opt_cl2_structure.write(structure,
                            work_path+rin.opt_cl2_cif,
                            symprec=rin.symtoleI,
                            title='ID_{0:d}'.format(next_id))
