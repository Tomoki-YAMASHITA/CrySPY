from logging import getLogger

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from ...util.struc_util import sort_by_atype, check_distance, get_nat, remove_zero
#from ...util.struc_util import cal_g, sort_by_atype_mol, find_site


logger = getLogger('cryspy')


def gen_permutation(
        atype,
        mindist,
        struc_data,
        sp,
        n_perm,
        id_start=None,
        symprec=0.01,
        ntimes=1,
        maxcnt_ea=50,
        struc_mol_id=None,
        molecular=False,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_perm (int): number of structures to generate by permutation
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    ntimes (int): number of swaps
    maxcnt_ea (int): maximum number of trial in permutation

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'permutation'}
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

    # ---------- generate structures by permutation
    while struc_cnt < n_perm:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ check nat for vc
        if len(parent_A.composition) == 1:
            logger.warning(f'Permutation: {pid_A} is composed of only single element')
            continue
        # ------ generate child
        if molecular:
            logger.error('Permutation for molecular is not implemented yet')
            # child, mol_id = pm.gen_child_mol(rin,
            #                                     self.sp.struc_data[pid_A],
            #                                     struc_mol_id[pid_A])
        else:
            child = gen_child(atype, mindist, parent_A, ntimes, maxcnt_ea)
        # ------ success
        if child is not None:
            children[cid] = child
            # if molecular:
            #     children_mol_id[cid] = mol_id
            parents[cid] = (pid_A, )    # tuple
            operation[cid] = 'permutation'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(child, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                    f' from {pid_A:>6} by permutation.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(atype, mindist, parent_A, ntimes=1, maxcnt_ea=50):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    ntimes (int): number of swaps
    maxcnt_ea (int): maximum number of trial in crossover

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize
    smatcher = StructureMatcher()    # instantiate StructureMatcher
    cnt = 0

    # ---------- ntimes permutation
    while True:
        child = parent_A.copy()    # keep original structure
        n = ntimes
        while n > 0:
            # ------ prepare index for each atom type
            indx_each_type = []
            for a in atype:
                indx_each_type.append(
                    [i for i, site in enumerate(child)
                        if site.species_string == a])
            # ------ choose two atom type
            non_empty_indices = [i for i, sublist in enumerate(indx_each_type) if sublist]
            type_choice = np.random.choice(non_empty_indices, 2, replace=False)
            # ------ choose index
            indx_choice = []
            for tc in type_choice:
                indx_choice.append(np.random.choice(indx_each_type[tc]))
            # ------ replace each other
            child.replace(indx_choice[0],
                                species=atype[type_choice[1]])
            child.replace(indx_choice[1],
                                species=atype[type_choice[0]])
            # ------ compare to original one
            if smatcher.fit(child, parent_A):
                n = ntimes    # back to the start
                continue
            else:
                n -= 1

        # ------ check distance
        success, mindist_ij, dist = check_distance(child, atype, mindist)
        if success:
            child = sort_by_atype(child, atype)
            return child
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in permutation: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Permutatin: could not satisfy min_dist' +
                        f' in {maxcnt_ea} times')
                logger.warning('Change parent')
                return None    # change parent
