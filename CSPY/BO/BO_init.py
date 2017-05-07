#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
import pandas as pd

from . import select_descriptor
from ..IO import pkl_data
from ..IO import read_input as rin


def initialize(stat, init_struc_data, rslt_data):
    print('\n#----------Initialize Bayesian optimization')
    with open('cspy.out', 'a') as fout:
        fout.write('\n#---------- Initilalize Bayesian optimization\n')

    #---------- initialize
    gen = 1
    id_done = np.array([], dtype=int)
    targets = np.array([], dtype=float)
    non_error_id = np.arange(len(init_struc_data))
    #------ rslt_data
    rslt_data['Gen'] = pd.Series(dtype=int)
    rslt_data = rslt_data[['Gen', 'Struc_ID', 'Spg_num', 'Spg_sym', 'Spg_num_opt',
                           'Spg_sym_opt', 'Energy', 'Magmom', 'Opt']]
    #- save
    pkl_data.save_rslt(rslt_data)


    #------ random select
    id_to_calc = random_select(len(init_struc_data), rin.interval)

    #---------- calc descriptor
    descriptors = select_descriptor.calc_X(init_struc_data)
    next_BO_id = len(init_struc_data)

    #------ save
    BO_id_data = (gen, next_BO_id, non_error_id, id_to_calc, id_done)
    pkl_data.save_BO_id(BO_id_data)
    BO_data = (descriptors, targets)
    pkl_data.save_BO_data(BO_data)

    #---------- status
    stat.set('status', 'generation', '{}'.format(gen))
    stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cspy.stat', 'w') as fstat:
        stat.write(fstat)

    #---------- out and log
    print('Generation: {}'.format(gen))
    print('selected_id: {}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cspy.out', 'a') as fout:
        fout.write('Generation: {}\n'.format(gen))
        fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))

    return rslt_data, BO_id_data, BO_data


def random_select(length, n):
    rnd_perm = np.random.permutation(xrange(length))
    selected_id = rnd_perm[0:n]
    return selected_id
