#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from ...IO import read_input as rin


def check_input_opt_cl2():
    #---------- prepare rin.jobfile, rin.opt_cl2_infile
    calc_inputs = [rin.jobfile, rin.opt_cl2_infile]

    #----- check required files
    for f in calc_inputs:
        if f == rin.opt_cl2_infile:
            finfiles = [rin.opt_cl2_infile + '_{}'.format(i) for i in range(1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    raise IOError('Could not find ./calc_in/'+ff)
        else:
            if not os.path.isfile('./calc_in/'+f):
                raise IOError('Could not find ./calc_in/'+f)


def clean_calc_files_opt_cl2(work_path):
    #---------- clean input files
    opt_cl2_files = [rin.opt_cl2_cif]
    for f in opt_cl2_files:
        if os.path.isfile(work_path+f):
            os.remove(work_path+f)

    #---------- clear stat file
    os.remove(work_path+'stat_job')
