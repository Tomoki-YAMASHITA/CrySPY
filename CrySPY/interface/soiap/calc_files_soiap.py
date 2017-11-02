#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from ...IO import read_input as rin


def check_input_soiap():
    # ---------- prepare rin.jobfile, rin.soiap_infile
    calc_inputs = [rin.jobfile, rin.soiap_infile]

    # ----- check required files
    for f in calc_inputs:
        if f == rin.soiap_infile:
            finfiles = [rin.soiap_infile + '_{}'.format(i) for i in range(1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    raise IOError('Could not find ./calc_in/'+ff)
        else:
            if not os.path.isfile('./calc_in/'+f):
                raise IOError('Could not find ./calc_in/'+f)


def clean_calc_files_soiap(work_path):
    # ---------- clean input files
    soiap_files = [rin.soiap_cif]
    for f in soiap_files:
        if os.path.isfile(work_path+f):
            os.remove(work_path+f)

    # ---------- clear stat file
    os.remove(work_path+'stat_job')
