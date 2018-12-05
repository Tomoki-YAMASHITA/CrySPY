#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np

from ..BO import combo_cryspy
from ..IO import pkl_data


def next_gen(stat, bo_id_data, bo_data):
    # ---------- out and log
    with open('cryspy.out', 'a') as fout:
        fout.write('# ---------- Bayesian optimization\n')
    print('# ---------- Bayesian optimization')

    # ---------- bo_id_data and bo_data
    gen, non_error_id, id_to_calc, id_done = bo_id_data
    descriptors, targets = bo_data

    # ---------- id_done --> sact
    sact = np.array([], dtype=int)
    for i in id_done:
        tindx = np.where(non_error_id == i)[0][0]
        sact = np.r_[sact, np.array([tindx])]

    # ---------- Bayesian optimization
    actions = combo_cryspy.bayes_opt(sact, descriptors, targets)

    # ---------- actions --> id_to_calc
    for i in actions:
        id_to_calc = np.r_[id_to_calc, non_error_id[i]]

    # ---------- gen
    gen += 1

    # ---------- save
    bo_id_data = (gen, non_error_id, id_to_calc, id_done)
    pkl_data.save_bo_id(bo_id_data)

    # ---------- status
    stat.set('status', 'generation', '{}'.format(gen))
    stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- out and log
    print('# ---------- Generation: {}'.format(gen))
    print('selected_id: {}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.out', 'a') as fout:
        fout.write('# ---------- Generation: {}\n'.format(gen))
        fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))
