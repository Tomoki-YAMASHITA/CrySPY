#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser
import random

import pandas as pd

from .select_descriptor import select_descriptor
from ..IO import pkl_data
from ..IO import read_input as rin


def initialize(stat, init_struc_data, rslt_data):
    # ---------- log
    print('\n# ---------- Selection: 1')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Selection 1\n')

    # ---------- check init_struc_data
    if None in init_struc_data.values():
        raise ValueError('init_struc_data includes None')

    # ---------- initialize
    n_selection = 1

    # ---------- rslt_data, add and sort
    rslt_data['Select'] = pd.Series(dtype=int)
    rslt_data = rslt_data[['Select', 'Struc_ID', 'Spg_num', 'Spg_sym', 'Spg_num_opt',
                           'Spg_sym_opt', 'E_eV_atom', 'Magmom', 'Opt']]
    pkl_data.save_rslt(rslt_data)

    # ---------- random select
    all_id = [i for i in range(len(init_struc_data))]
    if rin.manual_select_bo is not None:
        # ------ manual select bo
        print('Manual select: {}'.format(' '.join(str(i) for i in rin.manual_select_bo)))
        with open('cryspy.out', 'a') as fout:
            fout.write('Manual select: {}\n'.format(' '.join(str(i) for i in rin.manual_select_bo)))
        nselect = rin.interval - len(rin.manual_select_bo)
        id_to_calc = rin.manual_select_bo
        if 0 < nselect:
            diff_id = list(set(all_id) - set(id_to_calc))
            id_to_calc.extend(random.sample(diff_id, nselect))
        # -- delete the value for manual_select_bo in cryspy.in
        config = ConfigParser.ConfigParser()
        config.read('cryspy.in')
        config.set('BO', 'manual_select_bo', '')
        with open('cryspy.in', 'w') as f:
            config.write(f)
    else:
        id_to_calc = random.sample(all_id, rin.interval)

    # ---------- calc descriptor
    init_dscrpt_data = select_descriptor(init_struc_data)
    opt_dscrpt_data = {}  # initialize in dict

    # ---------- save for BO
    bo_id_data = (n_selection, id_to_calc)
    pkl_data.save_bo_id(bo_id_data)
    bo_data = (init_dscrpt_data, opt_dscrpt_data)
    pkl_data.save_bo_data(bo_data)

    # ---------- status
    stat.set('status', 'selection', '{}'.format(n_selection))
    stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- out and log
    print('selected_id: {}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.out', 'a') as fout:
        fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))
