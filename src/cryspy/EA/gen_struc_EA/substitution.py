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
        nat_data,
        ll_nat,
        ul_nat,
        id_start=None,
        symprec=0.01,
        maxcnt_ea=50,
        target='random',
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_subs (int): number of structures to generate by substitution
    nat_data (dict): {id: nat}
    ll_nat (tuple): lower limit of nat, e.g. (1, 1)
    ul_nat (tuple): upper limit of nat, e.g. (8, 8)
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    maxcnt_ea (int): maximum number of trial in substitution
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

    # ---------- generate structures by substitution
    while struc_cnt < n_subs:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ check nat limit
        atype_avail_elim = []
        atype_avail_add = []
        for i, at in enumerate(atype):
            if nat_data[pid_A][i] > ll_nat[i]:
                atype_avail_elim.append(at)
            if nat_data[pid_A][i] < ul_nat[i]:
                atype_avail_add.append(at)
        if len(atype_avail_elim) == 1:
            # if atype_avail_elim is ['Na'], 'Na' should be removed from atype_avail_add
            at = atype_avail_elim[0]
            if at in atype_avail_add:
                atype_avail_add.remove(at)
        if len(atype_avail_add) == 1:
            # if atype_avail_add is ['Na'], 'Na' should be removed from atype_avail_elim
            at = atype_avail_add[0]
            if at in atype_avail_elim:
                atype_avail_elim.remove(at)
        if len(atype_avail_add) == 0:
            logger.warning('Substitution: reached nat limit (ul_nat).')
            logger.warning('Change parent')
            continue
        if len(atype_avail_elim) == 0:
            logger.warning('Substitution: reached nat limit (ll_nat).')
            logger.warning('Change parent')
            continue
        child = gen_child(atype, mindist, parent_A, atype_avail_add, atype_avail_elim,
                          maxcnt_ea, target)
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


def gen_child(atype, mindist, parent_A, atype_avail_add, atype_avail_elim,
              maxcnt_ea=50, target='random'):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    atype_avail_add (list): available atom type for addition
    atype_avail_elim (list): available atom type for elimination
    maxcnt_ea (int): maximum number of trial in crossover
    target (str): only 'random' for now

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    cnt = 0

    # ---------- generate structure by substitution
    while True:
        child = parent_A.copy()    # keep original structure
        if target == 'random':
            # ------ random choice for atom type
            at = np.random.choice(atype_avail_elim)
            if at in atype_avail_add:
                atype_avail_add.remove(at)
            if len(atype_avail_add) == 0:
                logger.warning('Substitution: no atype_avail_add, retry.')
                cnt += 1
                if cnt >= maxcnt_ea:
                    logger.warning('Substitution: cnt >= maxcnt_ea.')
                    logger.warning('Change parent')
                    return None    # change parent
                continue
            else:
                # ------ elimination
                aindx = [i for i, site in
                            enumerate(child) if site.species_string == at]
                elim_indx = np.random.choice(aindx, 1)  # ", 1" to get array
                coords = child[elim_indx[0]].frac_coords    # coords of eliminated atom
                child.remove_sites(elim_indx)
                # ------ addition
                at = np.random.choice(atype_avail_add)
                child.append(species=at, coords=coords)
        # ---------- check mindist
        success, mindist_ij, dist = check_distance(child, atype, mindist)
        if success:
            child = sort_by_atype(child, atype)
            return child
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in addition: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Substitution: cnt >= maxcnt_ea.')
                logger.warning('Change parent')
                return None    # change parent
