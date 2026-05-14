from itertools import product
from logging import getLogger

import numpy as np

from ...util.struc_util import get_nat

# ---------- import later
#from ...util.charge_neutral import is_charge_neutral_nat


logger = getLogger('cryspy')


def gen_elimination(
        atype,
        struc_data,
        sp,
        n_elim,
        elim_dnat_map,
        id_start=None,
        symprec=0.01,
        target='random',
        rng=None,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_elim (int): number of structures to generate by elimination
    elim_dnat_map (dict): {id: list of delta nat combinations for elimination}
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    target (str): only 'random' for now
    rng (np.random.Generator): random number generator

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'strain'}
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    struc_cnt = 0
    children = {}
    #children_mol_id = {}
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

    # ---------- generate structures by elimination
    while struc_cnt < n_elim:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ delta nat conbinations
        dnat_comb = elim_dnat_map[pid_A]
        # ------ elim_element_list, e.g. ['Li', 'Li', 'O']
        if target == 'random':
            dnat = dnat_comb[rng.integers(len(dnat_comb))] 
            elim_element_list = [a for a, n in zip(atype, dnat) for _ in range(n)]
        # ------ generate child
        child = gen_child(parent_A, elim_element_list, rng=rng)
        # ------ success
        if child is not None:
            children[cid] = child
            parents[cid] = (pid_A, )
            operation[cid] = 'elimination'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(child, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {pid_A:>6} by elimination.'
                f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(parent_A, elim_element_list, rng=None):
    '''
    parent_A (Structure): pymatgen Structure object
    elim_element_list (list): list of atom types to eliminate, e.g. ['Li', 'Li', 'O']
    target (str): only 'random' for now
    rng (np.random.Generator): random number generator

    # ---------- return
    child (Structure): pymatgen Structure object
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- keep original structure
    child = parent_A.copy()    # keep original structure

    # ---------- elimination
    # ------ get indices of atoms to eliminate
    elim_indices = []
    used = set()
    for elem in elim_element_list:
        # まだ使われていないインデックスだけを候補に
        candidates = [i for i, site in enumerate(child) if site.species_string == elem and i not in used]
        idx = rng.choice(candidates)
        elim_indices.append(idx)
        used.add(idx)

    # ------ remove atom
    child.remove_sites(elim_indices)

    # ---------- return
    #            no need to check distance in elimination
    return child


def get_elim_dnat_comb(ll_nat, elim_max, parent_nat, cn_comb=None, charge=None):
    """
    Generate candidate delta nat (number of atoms) for elimination.

    This function returns all feasible atom-count decrements to be eliminated
    from `parent_nat` under the lower-bound and maximum-elimination constraints.

    Parameters
    ----------
    ll_nat : tuple[int]
        Per-element lower bounds of atom counts.
    elim_max : int
        Maximum total number of atoms to eliminate in one operation.
    parent_nat : tuple[int]
        nat of the parent structure.
    cn_comb (np.ndarray): charge neutral combinations
    charge (tuple): charge of each atom type

    Returns
    -------
    dnat_comb : list[tuple[int, ...]]
        Each candidate satisfies:
        - 0 < sum(dnat) <= elim_max
        - parent_nat - dnat >= ll_nat (element-wise)
        - sum(parent_nat - dnat) > 0
        - parent_nat - dnat is charge-neutral if charge is given and cn_comb is None

    Case 1: `cn_comb is None`
    All dnat candidates are locally enumerated under elim_max and ll_nat.

    Case 2: `cn_comb` is given
    Only `dnat` rows from `cn_comb` with sum(dnat) <= elim_max and
    parent_nat - dnat >= ll_nat are kept.

    Case 3: `cn_comb is None` and `charge` is given
    All dnat candidates are locally enumerated, and only candidates where
    parent_nat - dnat is charge-neutral are kept.
    """
    # ---------- charge-neutral checker
    if charge is not None and cn_comb is None:
        from ...util.charge_neutral import is_charge_neutral_nat

    # ---------- initialize
    dnat_comb = []

    # ---------- without precomputed charge-neutral combinations
    if cn_comb is None:
        max_del_per_element = [
            min(current - ll, elim_max)
            for current, ll in zip(parent_nat, ll_nat)
        ]

        for dnat in product(*[range(max_del + 1) for max_del in max_del_per_element]):
            new_nat = np.array(parent_nat) - np.array(dnat)
            if 0 < sum(dnat) <= elim_max and np.sum(new_nat) > 0 and np.all(new_nat >= ll_nat):
                # ------ charge-neutral check
                if charge is not None:
                    if not is_charge_neutral_nat(new_nat, charge):
                        continue

                dnat_comb.append(dnat)

    # ---------- with precomputed charge-neutral combinations
    else:
        mask = cn_comb.sum(axis=1) <= elim_max
        cn_comb_delta = cn_comb[mask].copy()    # delta combinations

        for dnat in cn_comb_delta:
            new_nat = np.array(parent_nat) - dnat
            if np.all(new_nat >= ll_nat) and np.sum(new_nat) > 0:
                dnat_comb.append(tuple(dnat))

    # ---------- return
    return dnat_comb
