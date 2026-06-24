from itertools import product
from logging import getLogger

import numpy as np

from ...util.struc_util import check_distance, sort_by_atype, get_nat


logger = getLogger('cryspy')


def gen_substitution(
        atype,
        mindist,
        struc_data,
        sp,
        n_subs,
        subs_comb_map,
        id_start=None,
        symprec=0.01,
        maxcnt_ea=50,
        target='random',
        rng=None
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_subs (int): number of structures to generate by substitution
    subs_comb_map (dict): {id: list of substitution patterns}
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    maxcnt_ea (int): maximum number of trial in substitution
    target (str): only 'random' for now
    rng (np.random.Generator): random number generator

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'substitution'}
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    struc_cnt = 0
    children = {}
    parents = {}
    operation = {}

    # ---------- id_offset
    if id_start is None:
        cid = max(struc_data.keys()) + 1
    else:
        if id_start < (max(struc_data.keys()) + 1):
            logger.error('id_start is already included in structure ID of the data')
        else:
            cid = id_start

    # ---------- generate structures by substitution
    while struc_cnt < n_subs:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ substitution combinations
        subs_comb = subs_comb_map[pid_A]
        # ------ subs_element_list, e.g. [('Li', 'Ca'), ('Ca', 'Cl')]
        if target == 'random':
            subs_element_list = subs_comb[int(rng.integers(len(subs_comb)))]
        # ------ generate child
        child = gen_child(atype, mindist, parent_A, subs_element_list, maxcnt_ea, target, rng=rng)
        # ------ success
        if child is not None:
            children[cid] = child
            parents[cid] = (pid_A, )
            operation[cid] = 'substitution'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(child, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {pid_A:>6} by substitution.'
                f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(atype, mindist, parent_A, subs_element_list, maxcnt_ea=50, target='random', rng=None):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    subs_element_list (list): list of tuples, e.g. [('Li', 'Ca'), ('Ca', 'Cl')]
    maxcnt_ea (int): maximum number of trial in crossover
    target (str): only 'random' for now
    rng (np.random.Generator): random number generator

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    cnt = 0

    # ---------- generate structure by substitution
    while True:
        child = parent_A.copy()    # keep original structure
        # ------ available indices for each atom type
        species_list = [str(site.specie) for site in child]
        available_indices = {elem: [i for i, sp in enumerate(species_list) if sp == elem] 
                            for elem in atype}
        # ------ select indices for substitution without duplicates
        subs_indices = []
        used_indices = set()
        for from_elem, to_elem in subs_element_list:
            available = [idx for idx in available_indices[from_elem] if idx not in used_indices]
            if not available:
                logger.error(f'Substitution invariant violated: no available indices for {from_elem}')
                raise RuntimeError(f'Substitution: no available indices for {from_elem}')
            selected_idx = rng.choice(available)
            subs_indices.append(selected_idx)
            used_indices.add(selected_idx)
        # ------ substitution
        for idx, (from_elem, to_elem) in zip(subs_indices, subs_element_list):
            child.replace(idx, to_elem)
        # ------ check mindist
        success, mindist_ij, dist = check_distance(child, atype, mindist)
        if success:
            child = sort_by_atype(child, atype)
            return child
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in substitution: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Substitution: cnt >= maxcnt_ea. Change parent.')
                return None    # change parent


def get_subs_comb(atype, ll_nat, ul_nat, subs_max, parent_nat, cn_set=None, charge=None):
    """
    Generate candidate substitution patterns.

    This function returns all feasible substitution patterns under the
    lower-bound, upper-bound, and maximum-substitution constraints.

    Parameters
    ----------
    atype : tuple[str]
        Atom types.
    ll_nat : tuple[int]
        Per-element lower bounds of atom counts.
    ul_nat : tuple[int]
        Per-element upper bounds of atom counts.
    subs_max : int
        Maximum total number of atoms to substitute in one operation.
    parent_nat : tuple[int]
        nat of the parent structure.
    cn_set : set[tuple[int, ...]]
        Set of charge-neutral nat tuples from precomputed cn_comb.
    charge : tuple
        Charge of each atom type.

    Returns
    -------
    subs_comb : list[list[tuple[str, str]]]
        List of substitution patterns.
        Each pattern is a list of (from_elem, to_elem) pairs.

    Notes
    -----
    If cn_set is given, charge neutrality is checked by membership in cn_set.
    If cn_set is None and charge is given, charge neutrality is checked from
    new_nat after substitution.
    """
    # ---------- charge-neutral checker
    if charge is not None and cn_set is None:
        from ...util.charge_neutral import is_charge_neutral_nat

    # ---------- initialize
    subs_comb = []

    # ---------- maximum number of subs for each element pair (no constraints)
    max_subs = {}
    for i, from_elem in enumerate(atype):
        for j, to_elem in enumerate(atype):
            if i != j:
                max_subs[(from_elem, to_elem)] = min(parent_nat[i], subs_max)

    # ---------- all combinations of substitution counts for each pair from 0 to max_subs
    pairs = list(max_subs.keys())
    max_counts = [max_subs[pair] for pair in pairs]

    # ---------- precompute element index and index pairs for speed
    idx = {elem: i for i, elem in enumerate(atype)}
    index_pairs = [(idx[from_elem], idx[to_elem]) for (from_elem, to_elem) in pairs]

    for counts in product(*[range(max_count + 1) for max_count in max_counts]):
        total_subs = sum(counts)
        if 0 < total_subs <= subs_max:
            # ------ new_nat after substitution
            new_nat = list(parent_nat)
            remain_nat = list(parent_nat)    # to track remaining atoms
            for (from_idx, to_idx), count in zip(index_pairs, counts):
                new_nat[from_idx] -= count
                remain_nat[from_idx] -= count
                if remain_nat[from_idx] < 0:
                    break
                new_nat[to_idx] += count

            # ------ skip if any element goes negative
            if any(n < 0 for n in remain_nat):
                continue

            # ------ check if new_nat is within limits
            if all(ll <= n <= ul for n, ll, ul in zip(new_nat, ll_nat, ul_nat)):
                # ------ charge-neutral check using precomputed cn_comb
                if cn_set is not None:
                    if tuple(new_nat) not in cn_set:
                        continue    # not charge-neutral --> skip this pattern

                # ------ charge-neutral check without precomputed cn_comb
                elif charge is not None:
                    if not is_charge_neutral_nat(new_nat, charge):
                        continue    # not charge-neutral --> skip this pattern

                # ------ OK: store this substitution pattern
                subs_list = []
                for pair, count in zip(pairs, counts):
                    subs_list.extend([pair] * count)
                if subs_list:
                    subs_comb.append(subs_list)

    # ---------- return
    return subs_comb
