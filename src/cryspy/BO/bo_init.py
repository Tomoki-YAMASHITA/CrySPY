'''
Initialize Bayesian optimization
'''

import configparser
from logging import getLogger
import random

import pandas as pd

from .select_descriptor import select_descriptor
from ..IO import io_stat, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def initialize(stat, init_struc_data, rslt_data):
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
    rslt_data = rslt_data[['Select', 'Spg_num',
                           'Spg_sym', 'Spg_num_opt',
                           'Spg_sym_opt', 'E_eV_atom', 'Magmom', 'Opt']]
    pkl_data.save_rslt(rslt_data)

    # ---------- random select
    all_id = [i for i in range(len(init_struc_data))]
    if rin.manual_select_bo:
        # ------ manual select bo
        x = ' '.join(str(i) for i in rin.manual_select_bo)
        logger.info(f'Manual select: {x}')
        nselect = rin.nselect_bo - len(rin.manual_select_bo)
        id_queueing = rin.manual_select_bo[:]    # shallow copy
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
    init_dscrpt_data = select_descriptor(init_struc_data)
    opt_dscrpt_data = {}  # initialize in dict

    # ---------- save for BO
    bo_id_data = (n_selection, id_queueing, id_running, id_select_hist)
    pkl_data.save_bo_id(bo_id_data)
    bo_data = (init_dscrpt_data, opt_dscrpt_data, bo_mean, bo_var, bo_score)
    pkl_data.save_bo_data(bo_data)

    # ---------- status
    io_stat.set_common(stat, 'selection', n_selection)
    io_stat.set_id(stat, 'selected_id', id_queueing)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- out and log
    x = ' '.join(str(a) for a in id_queueing)
    logger.info(f'selected_id: {x}')

