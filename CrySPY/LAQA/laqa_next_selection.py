#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

from ..IO import pkl_data
from ..IO import read_input as rin
from ..IO.out_results import out_laqa_id_hist


def next_selection(stat, laqa_id_data, laqa_data):
    # ---------- laqa_id_data and laqa_data
    (id_to_calc, id_select_hist, id_done) = laqa_id_data
    (tot_step_select, laqa_step, laqa_struc,
     laqa_energy, laqa_bias, laqa_score) = laqa_data

    # ---------- LAQA selection
    for k, v in sorted(laqa_score.items(), key=lambda x: -x[1][-1]):
        if v[-1] == -float('inf'):
            break
        else:
            id_to_calc.append(k)
            if len(id_to_calc) == rin.nselect:
                break

    # ---------- done LAQA
    if len(id_to_calc) == 0:
        with open('cryspy.out', 'a') as fout:
            fout.write('\nDone LAQA!\n')
        print('\nDone LAQA!')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- append id_select_hist and out
    id_select_hist.append(id_to_calc)
    out_laqa_id_hist(id_select_hist)

    # ---------- tot_step_select for next selection
    tot_step_select.append(0)

    # ---------- save
    laqa_id_data = (id_to_calc, id_select_hist, id_done)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc,
                 laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)

    # ---------- status
    stat.set('status', 'LAQA_selection', '{}'.format(len(id_select_hist)))
    if len(id_to_calc) > 30:
        stat.set('status', 'selected_id', '{} IDs'.format(len(id_to_calc)))
        stat.set('status', 'id_to_calc', '{} IDs'.format(len(id_to_calc)))
    else:
        stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
        stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- out and log
    print('\n# ---------- LAQA selection {}'.format(len(id_select_hist)))
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- LAQA selection {}\n'.format(len(id_select_hist)))
    if len(id_to_calc) > 30:
        print('selected_id: {} IDs'.format(len(id_to_calc)))
        with open('cryspy.out', 'a') as fout:
            fout.write('selected_id: {} IDs\n\n'.format(len(id_to_calc)))
    else:
        print('selected_id: {}\n'.format(' '.join(str(a) for a in id_to_calc)))
        with open('cryspy.out', 'a') as fout:
            fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))
