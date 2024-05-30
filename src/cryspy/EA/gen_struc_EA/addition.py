from logging import getLogger

import numpy as np

from ...util.struc_util import check_distance, sort_by_atype
#from .adj_comp import operation_atoms, convex_hull_check


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
        child = gen_child(atype, mindist, parent_A, atype_avail, maxcnt_ea, target)
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
            logger.info(f'Structure ID {cid:>6} was generated'
                f' from {pid_A:>6} by addition.'
                f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(atype, mindist, parent_A, atype_avail, maxcnt_ea=50, target='random'):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    atype_avail (list): available atom type for addition
    maxcnt_ea (int): maximum number of trial in crossover
    target (str): only 'random' for now

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    cnt = 0

    # ---------- generate child
    while True:
        child = parent_A.copy()    # keep original structure
        if target == 'random':
            # ------ random choice for atom type
            at = np.random.choice(atype_avail)
            # ------ add atom
            coords = np.random.rand(3)
            child.append(species=at, coords=coords)
        # ---------- not implemented yet
        # elif target in ['depop', 'overpop']:
        #     section = convex_hull_check(rin)
        #     child = operation_atoms(rin, 'addition', child, section)
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
                logger.warning('Addition: could not satisfy min_dist' +
                        f' in {maxcnt_ea} times')
                logger.warning('Change parent')
                return None    # change parent
