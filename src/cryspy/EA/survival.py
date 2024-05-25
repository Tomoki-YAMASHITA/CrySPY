from logging import getLogger

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher


logger = getLogger('cryspy')


def survival_fittest(fitness, struc_data, elite_struc=None, elite_fitness=None,
                     n_fittest=0, fit_reverse=False, emax_ea=None, emin_ea=None):
    '''
    # ---------- args
    fitness (dict): {ID: fitness, ...}
    struc_data (dict): {ID: structure data, ...}
    elite_struc (dict): {ID: elite structure data, ...}
    elite_fitness (dict): {ID: elite fitness, ...}
    n_fittest (int): number of fittest data which can survive
    fit_reverse (bool): if False, lower fitness is better
    emax_ea (float): maximum energy for cutoff
    emin_ea (float): minimum energy for cutoff

    # ---------- return
    ranking (list): [ID, ...]
    fit_with_elite (dict): {ID: fitness, ...} fitness + elite_fitness
    '''

    # ---------- initialize
    struc_with_elite = struc_data.copy()    # shallow copy to leave original data
    fit_with_elite = fitness.copy()          # shallow copy to leave original data
    ncheck = 5
    ranking = []
    smatcher = StructureMatcher()    # instantiate

    # ---------- elite
    if elite_struc is not None:
        struc_with_elite.update(elite_struc)
        fit_with_elite.update(elite_fitness)

    # ---------- None, np.nan, emax_ea, emin_ea
    for cid, value in fit_with_elite.items():
        if value is None or np.isnan(value):
            fit_with_elite[cid] = -np.inf if fit_reverse else np.inf
        # ------ emax_ea
        if emax_ea is not None:
            if value > emax_ea:
                fit_with_elite[cid] = -np.inf if fit_reverse else np.inf
                logger.info(f'Eliminate ID {cid}: {value} > emax_ea')
        # ------ emin_ea
        if emin_ea is not None:
            if value < emin_ea:
                fit_with_elite[cid] = -np.inf if fit_reverse else np.inf
                logger.info(f'Eliminate ID {cid}: {value} < emin_ea')

    # ---------- sort fit_with_elite
    sorted_fit_with_elite = sorted(fit_with_elite, key=fit_with_elite.get, reverse=fit_reverse)

    # ---------- ranking without duplication
    for cid in sorted_fit_with_elite:
        # ------ init dupl_flag
        dupl_flag = False
        # ------ for structure is None
        if struc_with_elite[cid] is None:
            continue    # next cid
        # ------ duplication check
        for tid in ranking[:-(ncheck+1):-1]:
            if smatcher.fit(struc_with_elite[cid], struc_with_elite[tid]):
                dupl_flag = True
                break
        # ------ register or skip
        if dupl_flag:
            continue    # next cid
        else:
            ranking.append(cid)
            # -- n_fittest
            if n_fittest > 0:
                if len(ranking) == n_fittest:
                    break

    # ---------- return
    return ranking, fit_with_elite, struc_with_elite