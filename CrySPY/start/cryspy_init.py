#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser
import os

import numpy as np
import pandas as pd

from . import utility
from ..gen_struc.random import rndgen
from ..IO import pkl_data
from ..IO import read_input as rin


def initialize():
    # ---------- start
    print(utility.get_date())
    print(utility.get_version())
    print('Start cryspy.py\n')
    with open('cryspy.out', 'w') as fout:
        fout.write(utility.get_date() + '\n')
        fout.write(utility.get_version() + '\n')
        fout.write('Start cryspy.py\n\n')

    # ---------- initialize stat
    stat = ConfigParser.ConfigParser()
    stat.add_section('input')
    stat.add_section('status')

    # ---------- read input
    print('Read input file, cryspy.in')
    rin.readin()          # read input data, cryspy,in
    rin.writeout()        # write input data in output file, cryspy.out
    rin.save_stat(stat)   # save input variables in cryspy.stat

    # ---------- make data directory
    if not os.path.isdir('data/pkl_data'):
        print('Make data directory')
        os.makedirs('data/pkl_data')

    # ---------- generate initial structures
    init_pos_path = utility.get_init_pos_path()
    if not rin.load_struc_flag:
        # ------ from scratch
        print('\n# --------- Generate initial structures')
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Generate initial structures\n')
        if rin.spgnum == 0:
            init_struc_data = rndgen.rndgen_wo_spg(
                                  rin.tot_struc, rin.natot, rin.atype, rin.nat, 0,
                                  rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                                  rin.maxcnt, rin.symtoleI, init_pos_path)
        else:
            fwpath = utility.check_fwpath()
            init_struc_data = rndgen.rndgen_spg(
                                  rin.tot_struc, rin.natot, rin.atype, rin.nat, rin.spgnum,
                                  0, rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                                  rin.maxcnt, rin.symtoleI, init_pos_path, fwpath)
        with open('cryspy.out', 'a') as fout:
            fout.write('Generated structures up to ID {}\n\n'.format(len(init_struc_data)-1))
        # ------ save
        pkl_data.save_init_struc(init_struc_data)

    else:
        # ------ load initial structure
        print('\n# --------- Load initial structure data')
        print('Load ./data/pkl_data/init_struc_data.pkl\n')
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Load initial structure data\n')
            fout.write('Load ./data/pkl_data/init_struc_data.pkl\n\n')
        init_struc_data = pkl_data.load_init_struc()
        # -- check
        if not rin.tot_struc == len(init_struc_data):
            raise ValueError('rin.tot_struc = {0}, len(init_struc_data) = {1}'.format(
                             rin.tot_struc, len(init_struc_data)))

    # ---------- initialize opt_struc_data
    opt_struc_data = []
    pkl_data.save_opt_struc(opt_struc_data)

    # ---------- initialize rslt_data
    rslt_data = pd.DataFrame(columns=['Struc_ID', 'Spg_num', 'Spg_sym',
                                      'Spg_num_opt', 'Spg_sym_opt',
                                      'Energy', 'Magmom', 'Opt'])
    rslt_data[['Struc_ID', 'Spg_num', 'Spg_num_opt']] = rslt_data[
                                   ['Struc_ID', 'Spg_num', 'Spg_num_opt']].astype(int)
    pkl_data.save_rslt(rslt_data)

    # ---------- return
    return stat, init_struc_data, opt_struc_data, rslt_data


def RS_init(stat):
    next_id = 0
    stat.set('status', 'next_id', '{}'.format(next_id))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)
    id_done = np.array([], dtype=int)
    RS_id_data = (next_id, id_done)
    pkl_data.save_RS_id(RS_id_data)

    return RS_id_data


def kpt_init():
    kpt_data = {}
    pkl_data.save_kpt(kpt_data)
    return kpt_data
