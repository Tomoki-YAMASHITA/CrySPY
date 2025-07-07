from itertools import product
from logging import getLogger
import random

import numpy as np

from ...util.struc_util import check_distance, sort_by_atype, get_nat


logger = getLogger('cryspy')


def gen_addition(
        atype,
        mindist,
        struc_data,
        sp,
        n_add,
        add_max,
        nat_data,
        ul_nat,
        id_start=None,
        symprec=0.01,
        maxcnt_ea=50,
        target='random',
        cn_comb_delta=None,
    ):
    '''
    # ---------- args
    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_add (int): number of structures to generate by addition
    add_max (int): maximum number of atoms to add
    nat_data (dict): {id: nat}
    ul_nat (tuple): upper limit of nat, e.g. (8, 8)
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    maxcnt_ea (int): maximum number of trial in addition
    target (str): target for addition, only 'random' for now
    cn_comb_delta (int): charge neutral combinations with a total number of atoms <= cn_nmax

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'addition'}
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

    # ---------- generate structures by addition
    while struc_cnt < n_add:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ delta nat conbinations
        dnat_comb =_get_dnat_comb(ul_nat, add_max, nat_data[pid_A], cn_comb_delta)
        if len(dnat_comb) == 0:
            logger.warning('Addition: no combinations found. Change parent')
            continue
        # ------ add_element_list, e.g. ['Li', 'Li', 'O']
        if target == 'random':
            dnat = random.choice(dnat_comb)
            add_element_list = [a for a, n in zip(atype, dnat) for _ in range(n)]
        # ------ generate child
        child = gen_child(
            atype,
            mindist,
            parent_A,
            add_element_list,
            maxcnt_ea,
        )
        # ------ success
        if child is not None:
            children[cid] = child
            parents[cid] = (pid_A, )
            operation[cid] = 'addition'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(child, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {pid_A:>6} by addition.'
                f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def _get_dnat_comb(ul_nat, add_max, parent_nat, cn_comb_delta):
    dnat_comb = []
    if cn_comb_delta is None:
        max_add_per_element = [min(ul - current, add_max) for ul, current in zip(ul_nat, parent_nat)]
        for dnat in product(*[range(max_add + 1) for max_add in max_add_per_element]):
            if 0 < sum(dnat) <= add_max:
                dnat_comb.append(dnat)
    else:    # charge neutrality
        for dnat in cn_comb_delta:
            new_nat = np.array(parent_nat) + dnat
            if np.all(new_nat <= ul_nat):
                dnat_comb.append(tuple(dnat))

    # ---------- return
    return dnat_comb


def gen_child(
        atype,
        mindist,
        parent_A,
        add_element_list,
        maxcnt_ea=50,
    ):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    add_element_list (list): list of atom types to add, e.g. ['Li', 'Li', 'O']
    maxcnt_ea (int): maximum number of trial in crossover
    target (str): only 'random' for now

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    vol10per = False
    vol20per = False
    child = parent_A.copy()    # keep original structure

    # ---------- generate child
    for at in add_element_list:
        cnt = 0
        while True:
            coords = np.random.rand(3)
            child.append(species=at, coords=coords)
            # ------ check distance
            success, mindist_ij, dist = check_distance(child, atype, mindist)
            if success:
                break
                # child = sort_by_atype(child, atype)
                # return child
            else:
                type0 = atype[mindist_ij[0]]
                type1 = atype[mindist_ij[1]]
                logger.warning(f'mindist in addition: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                child.pop()    # remove last atom
                if cnt >= maxcnt_ea:
                    # ------ volume change
                    if not vol10per:
                        vol10per = True
                        cnt = 0
                        child.scale_lattice(parent_A.volume * 1.1)
                        logger.warning('Addition: increase volume by 10%')
                        continue
                    elif vol10per and not vol20per:
                        vol20per = True
                        cnt = 0
                        child.scale_lattice(parent_A.volume * 1.2)
                        logger.warning('Addition: increase volume by 20%')
                        continue
                    # ------ fail
                    logger.warning('Addition: could not satisfy min_dist' +
                            f' in {maxcnt_ea} times')
                    logger.warning('Change parent')
                    return None    # change parent

    # ---------- complete
    child = sort_by_atype(child, atype)
    return child