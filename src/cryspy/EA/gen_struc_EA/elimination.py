from itertools import product
from logging import getLogger
import random

import numpy as np

from ...util.struc_util import get_nat


logger = getLogger('cryspy')


def gen_elimination(
        atype,
        struc_data,
        sp,
        n_elim,
        elim_max,
        nat_data,
        ll_nat,
        id_start=None,
        symprec=0.01,
        target='random',
        cn_comb_delta=None,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_elim (int): number of structures to generate by elimination
    elim_max (int): maximum number of atoms to eliminate
    nat_data (dict): {id: nat}
    ll_nat (tuple): lower limit of nat, e.g. (1, 1)
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    target (str): only 'random' for now
    cn_comb_delta (int): charge neutral combinations with a total number of atoms <= cn_nmax

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'strain'}
    '''

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
        dnat_comb = _get_dnat_comb(ll_nat, elim_max, nat_data[pid_A], cn_comb_delta)
        if len(dnat_comb) == 0:
            logger.warning('Elimination: no combinations found. Change parent')
            continue
        # ------ elim_element_list, e.g. ['Li', 'Li', 'O']
        if target == 'random':
            dnat = random.choice(dnat_comb)
            elim_element_list = [a for a, n in zip(atype, dnat) for _ in range(n)]
        # ------ generate child
        child = gen_child(parent_A, elim_element_list)
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


def _get_dnat_comb(ll_nat, elim_max, parent_nat, cn_comb_delta):
    dnat_comb = []
    if cn_comb_delta is None:
        max_del_per_element = [min(current, elim_max) for current in parent_nat]
        for dnat in product(*[range(max_del + 1) for max_del in max_del_per_element]):
            new_nat = np.array(parent_nat) - np.array(dnat)
            if 0 < sum(dnat) <= elim_max and np.sum(new_nat) > 0:
                dnat_comb.append(dnat)
    else:    # charge neutrality
        for dnat in cn_comb_delta:
            new_nat = np.array(parent_nat) - dnat
            if np.all(new_nat >= ll_nat) and np.sum(new_nat) > 0:
                dnat_comb.append(tuple(dnat))

    # ---------- return
    return dnat_comb


def gen_child(parent_A, elim_element_list):
    '''
    parent_A (Structure): pymatgen Structure object
    elim_element_list (list): list of atom types to eliminate, e.g. ['Li', 'Li', 'O']
    target (str): only 'random' for now

    # ---------- return
    child (Structure): pymatgen Structure object
    '''

    # ---------- keep original structure
    child = parent_A.copy()    # keep original structure

    # ---------- elimination
    # ------ get indices of atoms to eliminate
    elim_indices = []
    used = set()
    for elem in elim_element_list:
        # まだ使われていないインデックスだけを候補に
        candidates = [i for i, site in enumerate(child) if site.species_string == elem and i not in used]
        idx = random.choice(candidates)
        elim_indices.append(idx)
        used.add(idx)

    # ------ remove atom
    child.remove_sites(elim_indices)

    # ---------- return
    #            no need to check distance in elimination
    return child
