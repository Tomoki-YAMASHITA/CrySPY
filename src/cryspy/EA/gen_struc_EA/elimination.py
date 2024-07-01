from logging import getLogger

import numpy as np

from ...util.struc_util import get_nat


logger = getLogger('cryspy')


def gen_elimination(
        atype,
        struc_data,
        sp,
        n_elim,
        nat_data,
        ll_nat,
        id_start=None,
        symprec=0.01,
        target='random',
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_elim (int): number of structures to generate by elimination
    nat_data (dict): {id: nat}
    ll_nat (tuple): lower limit of nat, e.g. (1, 1)
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    target (str): only 'random' for now

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
        # ------ check nat limit
        atype_avail = []
        for i, at in enumerate(atype):
            if nat_data[pid_A][i] > ll_nat[i]:
                atype_avail.append(at)
        if len(atype_avail) == 0:
            logger.warning('Elimination: reached nat limit (ll_nat). cannot add atoms')
            logger.warning('Change parent')
            continue
        child = gen_child(parent_A, atype_avail, target)
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


def gen_child(parent_A, atype_avail, target='random'):
    '''
    parent_A (Structure): pymatgen Structure object
    atype_avail (list): available atom type for elimination
    target (str): only 'random' for now

    # ---------- return
    child (Structure): pymatgen Structure object
    '''

    # ---------- keep original structure
    child = parent_A.copy()    # keep original structure

    # ---------- elimination
    if target == 'random':
        # ------ random choice for atom type
        at = np.random.choice(atype_avail)
        # ------ random choice for atom index
        aindx = [i for i, site in enumerate(child) if site.species_string == at]
        elim_indx = np.random.choice(aindx, 1)  # ", 1" to get array
        # ------ remove atom
        child.remove_sites(elim_indx)

    # ---------- return
    #            no need to check distance in elimination
    return child
