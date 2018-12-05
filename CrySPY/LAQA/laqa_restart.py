#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..IO import pkl_data
from ..IO import read_input as rin


def restart(stat, laqa_id_data, laqa_data, prev_nstruc):
    # ---------- load laqa data
    id_to_calc, id_select_hist, id_done = laqa_id_data
    tot_step_select, laqa_step, laqa_struc, laqa_energy, laqa_bias, laqa_score = laqa_data

    # ---------- append scores and id_to_calc
    for i in range(prev_nstruc, rin.tot_struc):
        laqa_step[i] = []
        laqa_struc[i] = []
        laqa_energy[i] = []
        laqa_bias[i] = []
        laqa_score[i] = [float('inf')]
        id_to_calc.append(i)
    print('Append scores and id_to_calc')
    with open('cryspy.out', 'a') as fout:
        fout.write('Append scores and id_to_calc\n')

    # ---------- status
    if len(id_to_calc) > 30:
        stat.set('status', 'id_to_calc', '{} IDs'.format(len(id_to_calc)))
    else:
        stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- save for LAQA
    laqa_id_data = (id_to_calc, id_select_hist, id_done)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc, laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)
