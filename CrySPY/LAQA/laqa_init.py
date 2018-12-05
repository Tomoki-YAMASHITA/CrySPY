#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from ..IO import pkl_data
from ..IO import read_input as rin


def initialize(stat, init_struc_data):
    print('\n# ---------- Initialize LAQA')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Initilalize LAQA\n')

    # ---------- initialize
    tot_step_select = [0]
    laqa_step = {}
    laqa_struc = {}
    laqa_energy = {}
    laqa_bias = {}
    laqa_score = {}
    for i in range(rin.tot_struc):
        laqa_step[i] = []
        laqa_struc[i] = []
        laqa_energy[i] = []
        laqa_bias[i] = []
        laqa_score[i] = [float('inf')]
    id_to_calc = [i for i in range(rin.tot_struc)]
    id_select_hist = []
    id_done = []

    # ---------- save for LAQA
    laqa_id_data = (id_to_calc, id_select_hist, id_done)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc, laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)

    # ---------- status
    stat.set('status', 'LAQA_selection', '0')
    stat.set('status', 'total step', '0')
    if len(id_to_calc) > 30:
        stat.set('status', 'selected_id', '{} IDs'.format(len(id_to_calc)))
        stat.set('status', 'id_to_calc', '{} IDs'.format(len(id_to_calc)))
    else:
        stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
        stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- out and log
    print('# ---------- LAQA selection 0')
    with open('cryspy.out', 'a') as fout:
        fout.write('# ---------- LAQA selection 0\n')
    if len(id_to_calc) > 30:
        print('selected_id: {} IDs'.format(len(id_to_calc)))
        with open('cryspy.out', 'a') as fout:
            fout.write('selected_id: {} IDs\n\n'.format(len(id_to_calc)))
    else:
        print('selected_id: {}'.format(' '.join(str(a) for a in id_to_calc)))
        with open('cryspy.out', 'a') as fout:
            fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))
