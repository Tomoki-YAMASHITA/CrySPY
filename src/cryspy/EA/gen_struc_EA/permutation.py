from logging import getLogger

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from ...util.struc_util import sort_by_atype, check_distance, cal_g, sort_by_atype_mol, find_site


logger = getLogger('cryspy')


def gen_permutation(atype, mindist, struc_data, sp, n_perm,
                  id_start=None, symprec=0.001,
                  ntimes=1, maxcnt_ea=50,
                  struc_mol_id=None, molecular=False):
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
            logger.info(f'Structure ID {cid:>6} was generated'
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
            type_choice = np.random.choice(len(atype), 2,
                                            replace=False)
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




# class Permutation:
#     '''
#     permutation

#     # ---------- args
#     ntimes (int): ntimes permutaion

#     maxcnt_ea (int): maximum number of trial in crossover
#     # ---------- instance methods
#     self.gen_child(struc)
#         if success, return self.child
#         if fail, return None
#     '''

#     def __init__(self, mindist):
#         self.mindist = mindist

#     def gen_child(self, rin, struc):
#         '''
#         generate child structure

#         # ---------- return
#         (if success) self.child:
#         (if fail) None:
#         '''
#         # ---------- keep original structure
#         self.child = struc.copy()
#         # ---------- instantiate StructureMatcher
#         smatcher = StructureMatcher()    # instantiate
#         # ---------- ntimes permutation
#         cnt = 0
#         while True:
#             n = rin.ntimes
#             while n > 0:
#                 # ------ prepare index for each atom type
#                 indx_each_type = []
#                 for a in rin.atype:
#                     indx_each_type.append(
#                         [i for i, site in enumerate(self.child)
#                          if site.species_string == a])
#                 # ------ choose two atom type
#                 type_choice = np.random.choice(len(rin.atype), 2,
#                                                replace=False)
#                 # ------ choose index
#                 indx_choice = []
#                 for tc in type_choice:
#                     indx_choice.append(np.random.choice(indx_each_type[tc]))
#                 # ------ replace each other
#                 self.child.replace(indx_choice[0],
#                                    species=rin.atype[type_choice[1]])
#                 self.child.replace(indx_choice[1],
#                                    species=rin.atype[type_choice[0]])
#                 # ------ compare to original one
#                 if smatcher.fit(self.child, struc):
#                     n = rin.ntimes    # back to the start
#                     continue
#                 else:
#                     n -= 1
#             # ------ check distance
#             success, mindist_ij, dist = check_distance(self.child,
#                                                        rin.atype,
#                                                        self.mindist)
#             if success:
#                 self.child = sort_by_atype(self.child, rin.atype)
#                 return self.child
#             else:
#                 type0 = rin.atype[mindist_ij[0]]
#                 type1 = rin.atype[mindist_ij[1]]
#                 logger.warning(f'mindist in permutation: {type0} - {type1}, {dist}. retry.')
#                 cnt += 1
#                 if cnt >= rin.maxcnt_ea:
#                     logger.warning('Permutatin: could not satisfy min_dist' +
#                           f' in {rin.maxcnt_ea} times')
#                     logger.warning('Change parent')
#                     self.child = None
#                     return None    # change parent

#     def gen_child_mol(self, rin, struc, mol_id):
#         '''
#         generate child structures for mol

#         # ---------- return
#         (if success) self.child:
#         (if fail) None:
#         '''
#         # ---------- keep original structure
#         self.child = struc.copy()
#         self.group_id = mol_id[1]
#         self.mol_id = mol_id[0]
#         self.true_dists = mol_id[2]
#         # ---------- instantiate StructureMatcher
#         smatcher = StructureMatcher()    # instantiate
#         # ---------- ntimes permutation
#         cnt = 0
#         while True:
#             n = rin.ntimes
#             while n > 0:
#                 cnt += 1
#                 mol_g = cal_g(self.child, self.mol_id, self.group_id, self.true_dists)
#                 mol_type = list(set(self.group_id))
#                 self.child_frac_coords, self.group_id, self.child_species, self.mol_id = find_site(
#                     self.child,
#                     self.mol_id,
#                     self.group_id,
#                     self.true_dists)
#                 mol_choice = np.random.choice(mol_type, 2, replace=False)
#                 while self.check_mol_diff(mol_choice, self.mol_id, self.group_id):
#                     mol_choice = np.random.choice(mol_type, 2, replace=False)
#                 remove_index = []
#                 tmp_group_id = [tmp_group_id for tmp_group_id in self.group_id if tmp_group_id != mol_choice[0]]
#                 tmp_group_id = [tmp_group_id for tmp_group_id in tmp_group_id if tmp_group_id != mol_choice[1]]
#                 tmp_mol_id = []
#                 choiced_mol_id = []
#                 for mc in mol_choice:
#                     for i, gid in enumerate(self.group_id):
#                         if mc == gid:
#                             choiced_mol_id.append(self.mol_id[i])
#                             break
#                 for i, gi in enumerate(self.group_id):
#                     if mol_choice[0] != gi and mol_choice[1] != gi:
#                         tmp_mol_id.append(self.mol_id[i])
#                 self.child = Structure(self.child.lattice, self.child_species, self.child_frac_coords)
#                 for i, mid in enumerate(self.group_id):
#                     if mid == mol_choice[0]:
#                         coords = self.child_frac_coords[i] - mol_g[mol_choice[0]] + mol_g[mol_choice[1]]
#                         self.child.append(self.child_species[i], coords)
#                         remove_index.append(i)
#                         tmp_group_id.append(mol_choice[0])
#                         tmp_mol_id.append(choiced_mol_id[0])
#                     elif mid == mol_choice[1]:
#                         coords = self.child_frac_coords[i] - mol_g[mol_choice[1]] + mol_g[mol_choice[0]]
#                         self.child.append(self.child_species[i], coords)
#                         tmp_group_id.append(mol_choice[1])
#                         tmp_mol_id.append(choiced_mol_id[1])
#                         remove_index.append(i)
#                 self.group_id = tmp_group_id
#                 self.mol_id = tmp_mol_id
#                 self.child.remove_sites(remove_index)
#                 # ------ compare to original one
#                 if smatcher.fit(self.child, struc):
#                     n = rin.ntimes    # back to the start
#                     continue
#                 else:
#                     n -= 1
#             # ------ check distance
#             success, mindist_ij, dist = check_distance(self.child,
#                                                        rin.atype,
#                                                        self.mindist)
#             if success:
#                 self.child, self.mol_id, self.group_id = sort_by_atype_mol(self.child,
#                                                                            rin.atype,
#                                                                            self.mol_id,
#                                                                            self.group_id)
#                 return self.child, [self.mol_id, self.group_id, self.true_dists]
#             else:
#                 type0 = rin.atype[mindist_ij[0]]
#                 type1 = rin.atype[mindist_ij[1]]
#                 logger.warning(f'mindist in permutation: {type0} - {type1}, {dist}. retry.')
#                 cnt += 1
#                 if cnt >= rin.maxcnt_ea:
#                     logger.warning('Permutatin: could not satisfy min_dist' +
#                           f' in {rin.maxcnt_ea} times')
#                     logger.warning('Change parent')
#                     self.child = None
#                     return None    # change parent

#     def check_mol_diff(self, mol_choice, mol_id, group_id):
#         # ---------- check_mol_choice[0] and the mol_id
#         for i, gid in enumerate(group_id):
#             if gid == mol_choice[0]:
#                 choice_0_mol = mol_id[i]
#                 # ------ check_mol_choice[1] and the mol_id
#                 for i, gid in enumerate(group_id):
#                     if gid == mol_choice[1]:
#                         choice_1_mol = mol_id[i]
#                         # -- compare both mol_id
#                         if choice_0_mol == choice_1_mol:
#                             return True
#         return False
