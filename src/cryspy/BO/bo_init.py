import configparser
from logging import getLogger
import random

import pandas as pd

from .select_descriptor import select_descriptor
from ..IO import io_stat, pkl_data


logger = getLogger('cryspy')


def initialize(rin, init_struc_data, rslt_data):
    # ---------- log
    logger.info('# ---------- Selection: 1')

    # ---------- check init_struc_data
    if None in init_struc_data.values():
        raise ValueError('init_struc_data includes None')

    # ---------- initialize
    n_selection = 1
    id_running = []
    id_select_hist = []
    bo_mean = {}
    bo_var = {}
    bo_score = {}

    # ---------- rslt_data, add and sort
    rslt_data['Select'] = pd.Series(dtype=int)
    rslt_data = rslt_data[[
                    'Select',
                    'Spg_num',
                    'Spg_sym',
                    'Spg_num_opt',
                    'Spg_sym_opt',
                    'E_eV_atom',
                    'Magmom',
                    'Opt',
                ]]
    pkl_data.save_rslt(rslt_data)

    # ---------- random select
    all_id = [i for i in range(len(init_struc_data))]
    if rin.manual_select_bo:
        # ------ manual select bo
        x = ' '.join(str(i) for i in rin.manual_select_bo)
        logger.info(f'Manual select: {x}')
        nselect = rin.nselect_bo - len(rin.manual_select_bo)
        id_queueing = list(rin.manual_select_bo[:])    # shallow copy
        if 0 < nselect:
            diff_id = list(set(all_id) - set(id_queueing))
            id_queueing.extend(random.sample(diff_id, nselect))
        # ------ delete the value for manual_select_bo in cryspy.in
        config = configparser.ConfigParser()
        config.read('cryspy.in')
        config.set('BO', 'manual_select_bo', '')
        with open('cryspy.in', 'w') as f:
            config.write(f)
    else:
        id_queueing = random.sample(all_id, rin.nselect_bo)

    # ---------- id_select_hist
    id_select_hist.append(id_queueing[:])    # append shallow copy

    # ---------- calc descriptor
    init_dscrpt_data = select_descriptor(rin, init_struc_data)
    opt_dscrpt_data = {}  # initialize in dict

    # ---------- save for BO
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)
    pkl_data.save_n_selection(n_selection)
    pkl_data.save_id_select_hist(id_select_hist)
    pkl_data.save_init_dscrpt_data(init_dscrpt_data)
    pkl_data.save_opt_dscrpt_data(opt_dscrpt_data)
    pkl_data.save_bo_mean(bo_mean)
    pkl_data.save_bo_var(bo_var)
    pkl_data.save_bo_score(bo_score)

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_common(stat, 'selection', n_selection)
    io_stat.set_id(stat, 'selected_id', id_queueing)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- out and log
    x = ' '.join(str(a) for a in id_queueing)
    logger.info(f'selected_id: {x}')

