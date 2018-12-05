#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import pandas as pd

from ..IO import out_results
from ..IO import pkl_data
from ..IO import read_input as rin


def initialize(stat, rslt_data):
    # ---------- log
    print('\n# ---------- Initialize evolutionary algorithm')
    print('# ------ Generation 1')
    print('{} structures by random\n'.format(rin.tot_struc))
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Initilalize evolutionary algorithm\n')
        fout.write('# ------ Generation 1\n')
        fout.write('{} structures by random\n\n'.format(rin.tot_struc))

    # ---------- initialize
    next_id = 0
    id_done = np.array([], dtype=int)
    gen = 1
    # ------ ea_info
    ea_info = pd.DataFrame(columns=['Gen', 'Population',
                                    'Crossover', 'Permutation', 'Strain',
                                    'Random', 'Elite',
                                    'crs_func', 'crs_lat', 'slct_func'])
    ea_info.iloc[:, 0:7] = ea_info.iloc[:, 0:7].astype(int)
    tmp_info = pd.Series([1, rin.tot_struc, 0, 0, 0, rin.tot_struc, 0,
                          rin.crs_func, rin.crs_lat, rin.slct_func], index=ea_info.columns)
    ea_info = ea_info.append(tmp_info, ignore_index=True)
    out_results.out_ea_info(ea_info)
    # ------ ea_origin
    ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID', 'Operation', 'Parent'])
    ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)
    for cID in range(rin.tot_struc):
        tmp_origin = pd.Series([1, cID, 'random', None], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ elite
    elite_struc = None
    elite_fitness = None
    # ------ rslt_data
    rslt_data['Gen'] = pd.Series(dtype=int)
    rslt_data = rslt_data[['Gen', 'Struc_ID', 'Spg_num', 'Spg_sym', 'Spg_num_opt',
                           'Spg_sym_opt', 'Energy', 'Magmom', 'Opt']]

    # ---------- save
    ea_id_data = (gen, next_id, id_done)
    pkl_data.save_ea_id(ea_id_data)
    ea_data = (elite_struc, elite_fitness, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)
    pkl_data.save_rslt(rslt_data)

    # ---------- status
    stat.set('status', 'generation', '{}'.format(gen))
    stat.set('status', 'next_id', '{}'.format(next_id))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)
