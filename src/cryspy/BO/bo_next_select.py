import configparser
from contextlib import redirect_stdout
from logging import getLogger
import os

import numpy as np
import physbo

from ..IO import io_stat, out_results, pkl_data


logger = getLogger('cryspy')


def next_select(rin, rslt_data, bo_id_data, bo_data, noprint=False):
    # ---------- log
    logger.info('# ------ Bayesian optimization')

    # ---------- bo_id_data and bo_data
    n_selection, id_running, id_queueing, id_select_hist = bo_id_data
    init_dscrpt_data, opt_dscrpt_data, bo_mean, bo_var, bo_score = bo_data

    # ---------- n_selection
    n_selection += 1

    # ---------- manual select
    if rin.manual_select_bo is not None:
        manual_list = list(rin.manual_select_bo)
        x = ' '.join(str(i) for i in manual_list)
        logger.info(f'Manual select: {x}')
        # ------ check already selected
        tmp_list = [i for i in manual_list if i in opt_dscrpt_data]
        for j in tmp_list:
            manual_list.remove(j)
            logger.info(f'ID {j} was already selected. Remove it')
        # ------ number of structure to be selected
        nselect = rin.nselect_bo - len(manual_list)
        id_queueing = manual_list
        # ------ delete the value for manual_select_bo in cryspy.in
        config = configparser.ConfigParser()
        config.read('cryspy.in')
        config.set('BO', 'manual_select_bo', '')    # clear
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
                if rin.manual_select_bo is not None:
                    if i in rin.manual_select_bo:
                        # remove manualy selected IDs not to select them
                        non_error_id.remove(i)
                    else:
                        descriptors.append(init_dscrpt_data[i])
                else:
                    descriptors.append(init_dscrpt_data[i])
        # ------ targets
        targets = [rslt_data.loc[cid]['E_eV_atom'] for cid in done_id]
        # ------ list --> numpy
        s_act = np.array(s_act, dtype=int)
        descriptors = np.array(descriptors)
        targets = np.array(targets, dtype=float)
        # ------ Bayesian optimization
        actions, cryspy_mean, cryspy_var, cryspy_score = _bayes_opt(
            s_act, descriptors, targets, nselect,
            rin.score, rin.cdev, rin.num_rand_basis, noprint)
        # t_act, t_mean, t_var, t_score = _bayes_opt2(
        #     s_act, descriptors, targets, nselect,
        #     rin.score, rin.cdev, rin.num_rand_basis, noprint)
        # logger.debug(f't_act: {t_act}')
        # logger.debug(f't_mean: {t_mean}')
        # logger.debug(f't_var: {t_var}')
        # logger.debug(f't_score: {t_score}')
        # ------ actions --> id_queueing
        for i in actions:
            id_queueing.append(non_error_id[i])
        # ------ bo_mean, bo_var, bo_score
        remaining_id = list(set(non_error_id) - set(done_id))
        bo_mean[n_selection] = dict(zip(remaining_id, cryspy_mean))
        bo_var[n_selection] = dict(zip(remaining_id, cryspy_var))
        bo_score[n_selection] = dict(zip(remaining_id, cryspy_score))
        # ------ out
        out_results.out_bo_common('bo_mean', bo_mean, rin.tot_struc)
        out_results.out_bo_common('bo_var', bo_var, rin.tot_struc)
        out_results.out_bo_common('bo_score', bo_score, rin.tot_struc)
        out_results.out_bo_status(bo_mean, bo_var, bo_score,
                                  n_selection)
        # ------ save bo_data
        pkl_data.save_init_dscrpt_data(init_dscrpt_data)
        pkl_data.save_opt_dscrpt_data(opt_dscrpt_data)
        pkl_data.save_bo_mean(bo_mean)
        pkl_data.save_bo_var(bo_var)
        pkl_data.save_bo_score(bo_score)

    # ---------- id_select_hist
    id_select_hist.append(id_queueing[:])    # append shallow copy
    out_results.out_bo_id_hist(id_select_hist)

    # ---------- save
    pkl_data.save_id_queueing(id_queueing)
    # pkl_data.save_id_running(id_running)    # not used here
    pkl_data.save_n_selection(n_selection)
    pkl_data.save_id_select_hist(id_select_hist)

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_common(stat, 'selection', n_selection)
    io_stat.set_id(stat, 'selected_id', id_queueing)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- out and log
    logger.info(f'# ---------- Selection: {n_selection}')
    x = ' '.join(str(a) for a in id_queueing)
    logger.info(f'selected_id: {x}')


def _bayes_opt(
    s_act,
    descriptors,
    targets,
    nselect,
    score,
    cdev=0.001,
    num_rand_basis=0,
    noprint=False,
):
    '''
    # ---------- args
    s_act (ndarray): indexes of already calculated structures. NOT structure ID.
    descriptors (ndarray): descriptors of all structures, without error
    targets (ndarray): energies of already calculated structures
    nselect (int): number of structures to be selected
    score (str): 'PI', EI' or 'TS'
    cdev (float): cutoff value for standard deviation
    num_rand_basis (int): number of basis functions. if 0, gaussian process
    noprint (bool): if True, suppress print

    # ---------- return
    actions (ndarray): indexes of selected structures
    p_mean (ndarray): mean of predicted energies
    p_var (ndarray): variance of predicted energies
    p_score (ndarray): value of acquisition function
    '''

    # ---------- standardization
    dev = np.std(descriptors, axis=0)
    nonzero_indx = np.where(dev > cdev)[0]
    X = (descriptors[:, nonzero_indx] - np.mean(
        descriptors[:, nonzero_indx], axis=0)) / dev[nonzero_indx]

    # ---------- Declaring the policy by
    policy = physbo.search.discrete.policy(test_X=X, initial_data=[s_act, -targets])

    # ---------- Bayes search
    if noprint:
        with redirect_stdout(open(os.devnull, 'w')):
            actions = policy.bayes_search(
                                max_num_probes=1,
                                num_search_each_probe=nselect,
                                simulator=None,
                                score=score,
                                interval=0,
                                num_rand_basis = num_rand_basis,
                            )
    else:
        actions = policy.bayes_search(
                            max_num_probes=1,
                            num_search_each_probe=nselect,
                            simulator=None,
                            score=score,
                            interval=0,
                            num_rand_basis = num_rand_basis,
                        )

    # ---------- mean, var, score
    p_mean = policy.get_post_fmean(X)
    p_var = policy.get_post_fcov(X)
    p_score = policy.get_score(mode=score, xs=X)

    # ------ return
    return actions, p_mean, p_var, p_score