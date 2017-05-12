#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from ...IO import read_input as rin


def check_input_qe():
    #---------- prepare rin.jobfile, rin.infile, rin.pseudopots
    calc_inputs = [rin.jobfile, rin.infile] + rin.pseudopots

    #----- check required files
    for f in calc_inputs:
        if f == rin.infile:
            finfiles = [rin.infile + '_{}'.format(i) for i in range(1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/' + ff):
                    raise IOError('Could not find ./calc_in/' + ff)
        else:
            if not os.path.isfile('./calc_in/' + f):
                raise IOError('Could not find ./calc_in/' + f)


def clean_calc_files_qe(work_path):
    #---------- clean input files
    qe_files = ['tmp.cif']
    for f in qe_files:
        if os.path.isfile(work_path+f):
            os.remove(work_path+f)

    #---------- clear stat file
    os.remove(work_path+'stat_job')
