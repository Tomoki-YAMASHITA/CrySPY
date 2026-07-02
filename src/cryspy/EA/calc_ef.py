from logging import getLogger

import numpy as np


logger = getLogger('cryspy')


def calc_ef(energy, nat, end_point):
    '''
    Formation energy: eV/atom

    note:
    nat = [4, 3, 5]
    12 atoms: energy * 12 - 4 * end_point[0] - 3 * end_point[1] - 5 * end_point[2]
    1 atom: energy - 4/12 * end_point[0] - 3/12 * end_point[1] - 5/12 * end_point[2]
    --> energy - ratio[0] * end_point[0] - ratio[1] * end_point[1] - ratio[2] * end_point[2]
    '''

    # ---------- np.nan
    if np.isnan(energy):
        return np.nan

    # ---------- check
    if len(nat) != len(end_point):
        logger.error('len(nat) != len(end_point)')
        raise SystemExit(1)

    # ---------- calc formation energy
    ratio = tuple([x/sum(nat) for x in nat])
    ef = energy
    for x, y in zip(ratio, end_point):
        ef -= x * y

    # ---------- return
    return ef