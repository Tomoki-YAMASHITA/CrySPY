#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser

import numpy as np

from ..BO import combo_cryspy
from ..IO import pkl_data
from ..IO import read_input as rin


def next_select(stat, rslt_data, bo_id_data, bo_data):
    # ---------- out and log
    with open('cryspy.out', 'a') as fout:
        fout.write('# ------ Bayesian optimization\n')
    print('# ------ Bayesian optimization')

    # ---------- bo_id_data and bo_data
    n_selection, id_to_calc = bo_id_data
    init_dscrpt_data, opt_dscrpt_data = bo_data

    # ---------- manual select
    if rin.manual_select_bo is not None:
        print('Manual select: {}'.format(' '.join(str(i) for i in rin.manual_select_bo)))
        with open('cryspy.out', 'a') as fout:
            fout.write('Manual select: {}\n'.format(' '.join(str(i) for i in rin.manual_select_bo)))
        # ------ check already selected
        tmp_list = [i for i in rin.manual_select_bo if i in opt_dscrpt_data]
        for j in tmp_list:
            rin.manual_select_bo.remove(j)
            print('ID {} was already selected. Remove it'.format(j))
            with open('cryspy.out', 'a') as fout:
                fout.write('ID {} was already selected. Remove it\n'.format(j))
        # ------ number of structure to be selected
        nselect = rin.interval - len(rin.manual_select_bo)
        id_to_calc = rin.manual_select_bo
        # ------ delete the value for manual_select_bo in cryspy.in
        config = ConfigParser.ConfigParser()
        config.read('cryspy.in')
        config.set('BO', 'manual_select_bo', '')
        with open('cryspy.in', 'w') as f:
            config.write(f)
    else:
        nselect = rin.interval
    # ------ if number of remaining structures < nselect
    if len(init_dscrpt_data) - len(opt_dscrpt_data) < nselect:
        nselect = len(init_dscrpt_data) - len(opt_dscrpt_data)

    # ---------- selection
    if 0 < nselect:
        # ------ descriptors
        descriptors = []
        s_act = []
        done_id = []
        non_error_id = [i for i in range(len(init_dscrpt_data))]
        for i in range(len(init_dscrpt_data)):
            if i in opt_dscrpt_data:
                # -- already done
                if opt_dscrpt_data[i] is None:
                    non_error_id.remove(i)
                else:
                    s_act.append(len(descriptors))
                    done_id.append(i)
                    descriptors.append(opt_dscrpt_data[i])
            else:
                # -- not yet
                descriptors.append(init_dscrpt_data[i])
        # ------ targets
        targets = []
        for i in done_id:
            targets.append(rslt_data[rslt_data['Struc_ID'] == i]['E_eV_atom'].iat[0])
        # ------ list --> numpy
        s_act = np.array(s_act, dtype=int)
        descriptors = np.array(descriptors)
        targets = np.array(targets, dtype=float)
        # ------ Bayesian optimization
        actions = combo_cryspy.bayes_opt(s_act, descriptors, targets, nselect)
        # ------ actions --> id_to_calc
        for i in actions:
            id_to_calc.append(non_error_id[i])

    # ---------- n_selection
    n_selection += 1

    # ---------- save
    bo_id_data = (n_selection, id_to_calc)
    pkl_data.save_bo_id(bo_id_data)

    # ---------- status
    stat.set('status', 'selection', '{}'.format(n_selection))
    stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- out and log
    print('\n\n# ---------- Selection: {}'.format(n_selection))
    print('selected_id: {}'.format(' '.join(str(a) for a in id_to_calc)))
    with open('cryspy.out', 'a') as fout:
        fout.write('\n\n# ---------- Selection: {}\n'.format(n_selection))
        fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in id_to_calc)))
