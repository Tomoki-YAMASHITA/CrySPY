'''
Selection in Bayesian optimization
'''

import configparser
from contextlib import redirect_stdout
from logging import getLogger
import os

import numpy as np

from .combo_cryspy import Policy_cryspy
from ..IO import io_stat, out_results, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def next_select(stat, rslt_data, bo_id_data, bo_data, noprint=False):
    # ---------- log
    logger.info('# ------ Bayesian optimization')

    # ---------- bo_id_data and bo_data
    n_selection, id_running, id_queueing, id_select_hist = bo_id_data
    init_dscrpt_data, opt_dscrpt_data, bo_mean, bo_var, bo_score = bo_data

    # ---------- n_selection
    n_selection += 1

    # ---------- manual select
    if rin.manual_select_bo:
        x = ' '.join(str(i) for i in rin.manual_select_bo)
        logger.info(f'Manual select: {x}')
        # ------ check already selected
        tmp_list = [i for i in rin.manual_select_bo if i in opt_dscrpt_data]
        for j in tmp_list:
            rin.manual_select_bo.remove(j)
            logger.info(f'ID {j} was already selected. Remove it')
        # ------ number of structure to be selected
        nselect = rin.nselect_bo - len(rin.manual_select_bo)
        id_queueing = rin.manual_select_bo
        # ------ delete the value for manual_select_bo in cryspy.in
        config = configparser.ConfigParser()
        config.read('cryspy.in')
        config.set('BO', 'manual_select_bo', '')
        with open('cryspy.in', 'w') as f:
            config.write(f)
    else:
        nselect = rin.nselect_bo
    # ------ if number of remaining structures < nselect
    if len(init_dscrpt_data) - len(opt_dscrpt_data) < nselect:
        nselect = len(init_dscrpt_data) - len(opt_dscrpt_data)

    # ---------- selection
    if 0 < nselect:
        # ------ descriptors
        descriptors = []
        s_act = []
        done_id = []    # finished and non_error
        non_error_id = [i for i in range(len(init_dscrpt_data))]
        for i in range(len(init_dscrpt_data)):
            if i in opt_dscrpt_data:
                # -- already done
                if opt_dscrpt_data[i] is None:    # find error
                    non_error_id.remove(i)
                    continue
                x = rslt_data.loc[i]['E_eV_atom']
                if rin.emax_bo is not None:
                    if x > rin.emax_bo:
                        non_error_id.remove(i)
                        logger.info(f'Eliminate ID {i}: {x} > emax_bo')
                        continue
                if rin.emin_bo is not None:
                    if x < rin.emin_bo:
                        non_error_id.remove(i)
                        logger.info(f'Eliminate ID {i}: {x} < emin_bo')
                        continue
                s_act.append(len(descriptors))
                done_id.append(i)
                descriptors.append(opt_dscrpt_data[i])
            else:
                # -- not yet
                if i in rin.manual_select_bo:
                    # remove manualy selected IDs not to select them
                    non_error_id.remove(i)
                else:
                    descriptors.append(init_dscrpt_data[i])
        # ------ targets
        targets = [rslt_data.loc[cid]['E_eV_atom'] for cid in done_id]
        # ------ list --> numpy
        s_act = np.array(s_act, dtype=int)
        descriptors = np.array(descriptors)
        targets = np.array(targets, dtype=float)
        # ------ Bayesian optimization
        actions, cryspy_mean, cryspy_var, cryspy_score = bayes_opt(
            s_act, descriptors, targets, nselect, noprint)
        # ------ actions --> id_queueing
        for i in actions:
            id_queueing.append(non_error_id[i])
        # ------ bo_mean, bo_var, bo_score
        remaining_id = list(set(non_error_id) - set(done_id))
        bo_mean[n_selection] = dict(zip(remaining_id, cryspy_mean))
        bo_var[n_selection] = dict(zip(remaining_id, cryspy_var))
        if rin.score == 'TS':
            # -- in TS, nested list [[score1, socre2, ...]] is used
            bo_score[n_selection] = dict(zip(remaining_id, cryspy_score[0]))
        else:
            bo_score[n_selection] = dict(zip(remaining_id, cryspy_score))
        # ------ out
        out_results.out_bo_common('BO_mean', bo_mean, rin.tot_struc)
        out_results.out_bo_common('BO_var', bo_var, rin.tot_struc)
        out_results.out_bo_common('BO_score', bo_score, rin.tot_struc)
        out_results.out_bo_status(bo_mean, bo_var, bo_score,
                                  n_selection)
        # ------ save bo_data
        bo_data = (init_dscrpt_data, opt_dscrpt_data,
                   bo_mean, bo_var, bo_score)
        pkl_data.save_bo_data(bo_data)

    # ---------- id_select_hist
    id_select_hist.append(id_queueing[:])    # append shallow copy
    out_results.out_bo_id_hist(id_select_hist)

    # ---------- save
    bo_id_data = (n_selection, id_queueing, id_running, id_select_hist)
    pkl_data.save_bo_id(bo_id_data)

    # ---------- status
    io_stat.set_common(stat, 'selection', n_selection)
    io_stat.set_id(stat, 'selected_id', id_queueing)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- out and log
    logger.info(f'# ---------- Selection: {n_selection}')
    x = ' '.join(str(a) for a in id_queueing)
    logger.info(f'selected_id: {x}')

    # ---------- ext
    if rin.calc_code == 'ext':
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('out\n')


def bayes_opt(s_act, descriptors, targets, nselect, noprint=False):
    # ---------- start COMBO part
    # ------ standardization
    # X = combo.misc.centering(descriptors)
    dev = np.std(descriptors, axis=0)
    nonzero_indx = np.where(dev > rin.cdev)[0]
    X = (descriptors[:, nonzero_indx] - np.mean(
        descriptors[:, nonzero_indx], axis=0)) / dev[nonzero_indx]

    # ------ Declaring the policy by
    pc = Policy_cryspy(test_X=X)

    # ------ pick up data, already optimized
    actions = pc.specified_search(specified_actions=s_act, max_num_probes=1)

    # ------ write
    pc.write(actions, -targets)    # minus for minimum search

    # ------ log
    out_log(pc.history)

    # ------ Bayes_search
    if noprint:
        with redirect_stdout(open(os.devnull, 'w')):
            actions, cryspy_mean, cryspy_var, cryspy_score = pc.bayes_search_cryspy(
                max_num_probes=1,
                num_search_each_probe=nselect,
                score=rin.score,
                num_rand_basis=rin.num_rand_basis)
    else:
        actions, cryspy_mean, cryspy_var, cryspy_score = pc.bayes_search_cryspy(
            max_num_probes=1,
            num_search_each_probe=nselect,
            score=rin.score,
            num_rand_basis=rin.num_rand_basis)

    # ------ return
    return actions, cryspy_mean, cryspy_var, cryspy_score


# ---------- just for log
def out_log(history):
    n = history.total_num_search
    index = np.argmax(history.fx[0:n])

    logger.info(f'current best E = {-history.fx[index]}')
