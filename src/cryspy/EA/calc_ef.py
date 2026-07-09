from logging import getLogger

import numpy as np


logger = getLogger('cryspy')


def calc_ef(energy, nat, ref_energies):
    '''
    Formation energy: eV/atom

    note:
    nat = [4, 3, 5]
    12 atoms: energy * 12 - 4 * ref_energies[0] - 3 * ref_energies[1] - 5 * ref_energies[2]
    1 atom: energy - 4/12 * ref_energies[0] - 3/12 * ref_energies[1] - 5/12 * ref_energies[2]
    --> energy - ratio[0] * ref_energies[0] - ratio[1] * ref_energies[1] - ratio[2] * ref_energies[2]
    '''

    # ---------- np.nan
    if np.isnan(energy):
        return np.nan

    # ---------- check
    if len(nat) != len(ref_energies):
        logger.error('len(nat) != len(ref_energies)')
        raise SystemExit(1)

    # ---------- calc formation energy
    ratio = tuple([x/sum(nat) for x in nat])
    ef = energy
    for x, y in zip(ratio, ref_energies):
        ef -= x * y

    # ---------- return
    return ef