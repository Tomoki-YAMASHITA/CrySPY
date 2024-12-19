from logging import getLogger

import numpy as np

from ...util.struc_util import check_distance, sort_by_atype, get_nat
from ...util.struc_util import get_cn_comb_within_n

logger = getLogger('cryspy')


def gen_addition(
        atype,
        mindist,
        struc_data,
        sp,
        n_add,
        nat_data,
        ul_nat,
        id_start=None,
        symprec=0.01,
        maxcnt_ea=50,
        target='random',
        charge=None,
        cn_nmax=3,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_add (int): number of structures to generate by addition
    nat_data (dict): {id: nat}
    ul_nat (tuple): upper limit of nat, e.g. (8, 8)
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    maxcnt_ea (int): maximum number of trial in addition

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
        # ------ check nat limit
        atype_avail = []
        for i, at in enumerate(atype):
            if nat_data[pid_A][i] < ul_nat[i]:
                atype_avail.append(at)
        if len(atype_avail) == 0:
            logger.warning('Addition: reached nat limit (ul_nat). cannot add atoms')
            logger.warning('Change parent')
            continue
        # ------ charge neutrality
        if charge is not None:
            cn_comb = get_cn_comb_within_n(charge, cn_nmax)
        else:
            cn_comb = None
        # ------ generate child
        child = gen_child(
            atype,
            mindist,
            parent_A,
            atype_avail,
            maxcnt_ea,
            target,
            cn_comb,
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


def gen_child(
        atype,
        mindist,
        parent_A,
        atype_avail,
        maxcnt_ea=50,
        target='random',
        cn_comb=None,
    ):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    atype_avail (list): available atom type for addition
    maxcnt_ea (int): maximum number of trial in crossover
    target (str): only 'random' for now
    cn_comb (tuple): charge neutral combinations of atoms within n atoms

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    cnt = 0
    vol10per = False
    vol20per = False

    # ---------- generate child
    while True:
        child = parent_A.copy()    # keep original structure
        if vol10per and not vol20per:
            child.scale_lattice(parent_A.volume * 1.1)
        elif vol20per:
            child.scale_lattice(parent_A.volume * 1.2)
        if target == 'random':
            # ------ random choice for atom type
            at = np.random.choice(atype_avail)
            # ------ add atom
            coords = np.random.rand(3)
            child.append(species=at, coords=coords)
        # ------ check distance
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
                # ------ volume change
                if not vol10per:
                    vol10per = True
                    cnt = 0    # reset cnt
                    logger.warning('Addition: increase volume by 10%')
                    continue
                elif vol10per and not vol20per:
                    vol20per = True
                    cnt = 0
                    logger.warning('Addition: increase volume by 20%')
                    continue
                # ------ fail
                logger.warning('Addition: could not satisfy min_dist' +
                        f' in {maxcnt_ea} times')
                logger.warning('Change parent')
                return None    # change parent
