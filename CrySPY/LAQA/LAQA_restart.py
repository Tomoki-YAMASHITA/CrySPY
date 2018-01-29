#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..IO import pkl_data
from ..IO import read_input as rin


def restart(stat, LAQA_id_data, LAQA_data, prev_nstruc):
    # ---------- load LAQA data
    id_to_calc, id_select_hist, id_done = LAQA_id_data
    total_step, LAQA_step, LAQA_struc, LAQA_energy, LAQA_bias, LAQA_score = LAQA_data

    # ---------- append scores and id_to_calc
    for i in range(prev_nstruc, rin.tot_struc):
        LAQA_step[i] = []
        LAQA_struc[i] = []
        LAQA_energy[i] = []
        LAQA_bias[i] = []
        LAQA_score[i] = [float('inf')]
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
    LAQA_id_data = (id_to_calc, id_select_hist, id_done)
    pkl_data.save_LAQA_id(LAQA_id_data)
    LAQA_data = (total_step, LAQA_step, LAQA_struc, LAQA_energy, LAQA_bias, LAQA_score)
    pkl_data.save_LAQA_data(LAQA_data)

    # ---------- return
    return stat, LAQA_id_data, LAQA_data
