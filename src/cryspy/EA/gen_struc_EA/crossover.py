'''
Crossover class
'''

from collections import Counter
from logging import getLogger
import random

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import DummySpecie

from ...IO import read_input as rin
from ...util.struc_util import origin_shift, sort_by_atype, check_distance
from ...util.struc_util import find_site, cal_g, sort_by_atype_mol, get_nat


logger = getLogger('cryspy')

class Crossover:
    '''
    crossover

    # ---------- instance methods
    self.gen_child(struc_A, struc_B)
        if success, return self.child
        if fail, return None
    '''

    def __init__(self, mindist):
        # ---------- self
        self.mindist = mindist

        # ---------- weight of lattice
        if rin.crs_lat == 'equal':
            self.w_lat = np.array([1.0, 1.0])
        elif rin.crs_lat == 'random':
            self.w_lat = np.random.choice([0.0, 1.0], size=2, replace=False)
        else:
            logger.error('crs_lat must be equal or random')

    def gen_child(self, struc_A, struc_B):
        '''
        generate child structure

        # ---------- return
        (if success) self.child:
        (if fail) None:
        '''
        # ---------- initialize
        self.parent_A = origin_shift(struc_A)
        self.parent_B = origin_shift(struc_B)
        count = 0
        # ---------- lattice crossover
        self._lattice_crossover()
        # ---------- generate children
        while True:
            count += 1
            # ------ coordinate crossover
            self._one_point_crossover()
            self.child = Structure(lattice=self.lattice, species=self.species,
                                   coords=self.coords)
            # ------ check nat_diff
            self._get_nat_diff()    # get self._nat_diff
            if any([abs(n) > rin.nat_diff_tole for n in self._nat_diff]):
                if count > rin.maxcnt_ea:    # fail
                    self.child = None
                    return self.child
                continue    # slice again
            # ------ check mindist
            dist_list = check_distance(self.child, rin.atype,
                                       self.mindist, check_all=True)
            # ------ something smaller than mindist
            if dist_list:
                # -- remove atoms within mindist
                if any([n > 0 for n in self._nat_diff]):
                    self._remove_within_mindist()
                    if self.child is None:    # fail --> slice again
                        if count > rin.maxcnt_ea:
                            return None
                        continue
                else:    # nothing to remove, nat_diff = [0, 0]
                    if count > rin.maxcnt_ea:
                        return None
                    continue    # fail --> slice again
            # ------ recheck nat_diff
            self._get_nat_diff()
            # ------ nothing smaller than mindist
            # -- remove atoms near the border line
            if any([n > 0 for n in self._nat_diff]):
                self._remove_border_line()
            # -- add atoms near border line
            if any([n < 0 for n in self._nat_diff]):
                self._add_border_line()
            # -- success --> break while loop
            if self.child is not None:
                break
            # -- fail --> slice again
            else:
                if count > rin.maxcnt_ea:
                    return None
                continue
        # ---------- final check for nat
        self._get_nat_diff()
        if not all([n == 0 for n in self._nat_diff]):
            return None    # failure
        # ---------- sort by atype
        self.child = sort_by_atype(self.child, rin.atype)
        # ---------- return
        return self.child

    def _lattice_crossover(self):
        # ---------- component --> self.w_lat
        matrix = ((self.w_lat[0]*self.parent_A.lattice.matrix
                  + self.w_lat[1]*self.parent_B.lattice.matrix)
                  / self.w_lat.sum())
        mat_len = np.sqrt((matrix**2).sum(axis=1))
        # ---------- absolute value of vector
        lat_len = ((np.array(self.parent_A.lattice.abc)*self.w_lat[0]
                   + np.array(self.parent_B.lattice.abc)*self.w_lat[1])
                   / self.w_lat.sum())
        # ---------- correction of vector length
        lat_array = np.empty([3, 3])
        for i in range(3):
            lat_array[i] = matrix[i]*lat_len[i]/mat_len[i]
        # ---------- Lattice for pymatgen
        self.lattice = Lattice(lat_array)

    def _one_point_crossover(self):
        # ---------- slice point
        while True:
            self._slice_point = np.random.normal(loc=0.5, scale=0.1)
            if 0.3 <= self._slice_point <= 0.7:
                break
        self._axis = np.random.choice([0, 1, 2])
        # ---------- crossover
        species_A = []
        species_B = []
        coords_A = []
        coords_B = []
        for i in range(self.parent_A.num_sites):
            if self.parent_A.frac_coords[i, self._axis] <= self._slice_point:
                species_A.append(self.parent_A[i].species_string)
                coords_A.append(self.parent_A[i].frac_coords)
            else:
                species_B.append(self.parent_A[i].species_string)
                coords_B.append(self.parent_A[i].frac_coords)
            if self.parent_B.frac_coords[i, self._axis] >= self._slice_point:
                species_A.append(self.parent_B[i].species_string)
                coords_A.append(self.parent_B[i].frac_coords)
            else:
                species_B.append(self.parent_B[i].species_string)
                coords_B.append(self.parent_B[i].frac_coords)
        # ---------- adopt a structure with more atoms
        if len(species_A) > len(species_B):
            species = species_A
            coords = coords_A
        elif len(species_A) < len(species_B):
            species = species_B
            coords = coords_B
        else:
            if np.random.choice([0, 1]):
                species = species_A
                coords = coords_A
            else:
                species = species_B
                coords = coords_B
        # ---------- set instance variables
        self.species, self.coords = species, coords

    def _get_nat_diff(self):
        tmp_nat, _ = get_nat(self.child, rin.atype)
        self._nat_diff = [i - j for i, j in zip(tmp_nat, rin.nat)]
        
    def _remove_within_mindist(self):
        '''
        if success: self.child <-- child structure data
        if fail: self.child <-- None
        '''
        for itype in range(len(rin.atype)):
            while self._nat_diff[itype] > 0:
                # ---------- check dist
                dist_list = check_distance(self.child, rin.atype,
                                           self.mindist, check_all=True)
                if not dist_list:    # nothing within mindist
                    return
                # ---------- appearance frequency
                ij_within_dist = [isite[0] for isite in dist_list] + [
                    jsite[1] for jsite in dist_list]
                site_counter = Counter(ij_within_dist)
                # ---------- get index for removing
                rm_index = None
                # ---- site[0]: index, site[1]: count
                for site in site_counter.most_common():
                    if self.child[site[0]].species_string == rin.atype[itype]:
                        rm_index = site[0]
                        break    # break for loop
                # ---------- remove atom
                if rm_index is None:
                    self.child = None
                    return
                else:
                    self.child.remove_sites([rm_index])
                    self._nat_diff[itype] -= 1
        # ---------- final check
        dist_list = check_distance(self.child, rin.atype,
                                   self.mindist, check_all=True)
        if dist_list:    # still something within mindist
            self.child = None
            logger.warning('some atoms within mindist. retry.')

    def _remove_border_line(self):
        # ---------- rank atoms from border line
        coords_axis = self.child.frac_coords[:, self._axis]
        # ------ one point crossover: boundary --> 0.0, slice_point, 1.0
        near_sp = (self._slice_point/2.0 < coords_axis) & \
            (coords_axis < (self._slice_point + 1.0)/2.0)
        near_one = (self._slice_point + 1.0)/2.0 <= coords_axis
        # -- distance from nearest boundary
        coords_diff = np.where(near_sp,
                               abs(coords_axis - self._slice_point),
                               coords_axis)
        coords_diff = np.where(near_one, 1.0 - coords_diff, coords_diff)
        atom_border_indx = np.argsort(coords_diff)
        # ---------- remove list
        rm_list = []
        for itype, nrm in enumerate(self._nat_diff):
            rm_list.append([])
            if nrm > 0:
                for ab_indx in atom_border_indx:
                    if self.child[ab_indx].species_string == rin.atype[itype]:
                        rm_list[itype].append(ab_indx)
                    if len(rm_list[itype]) == nrm:
                        break
        # ---------- remove
        for each_type in rm_list:
            if each_type:
                self.child.remove_sites(each_type)

    def _add_border_line(self):
        for i in range(len(rin.atype)):
            # ---------- counter
            cnt = 0
            # ---------- add atoms
            while self._nat_diff[i] < 0:
                cnt += 1
                coords = np.random.rand(3)
                self._mean_choice()
                coords[self._axis] = np.random.normal(loc=self._mean,
                                                      scale=0.08)
                self.child.append(species=rin.atype[i], coords=coords)
                success, mindist_ij, dist = check_distance(self.child,
                                                           rin.atype,
                                                           self.mindist)
                if success:
                    cnt = 0    # reset
                    self._nat_diff[i] += 1
                else:
                    type0 = rin.atype[mindist_ij[0]]
                    type1 = rin.atype[mindist_ij[1]]
                    logger.warning(f'mindist in _add_border_line: {type0} - {type1}, {dist}. retry.')
                    self.child.pop()    # cancel
                # ------ fail
                if cnt == rin.maxcnt_ea:
                    self.child = None
                    return

    def _mean_choice(self):
        '''which boundary possesses more atoms'''
        n_zero = np.sum(np.abs(self.child.frac_coords[:, self._axis] - 0.0)
                        < 0.1)
        n_slice = np.sum(np.abs(self.child.frac_coords[:, self._axis]
                                - self._slice_point) < 0.1)
        if n_zero < n_slice:
            self._mean = 0.0
        elif n_zero > n_slice:
            self._mean = self._slice_point
        else:
            self._mean = np.random.choice([0.0, self._slice_point])

    def gen_child_mol(self, struc_A, struc_B, mol_id_A, mol_id_B):
        '''
        generate child structures for mol

        # ---------- return
        (if success) self.child:
        (if fail) None:
        '''
        # ---------- initialize
        self.parent_A = origin_shift(struc_A)
        self.parent_B = origin_shift(struc_B)
        self.mol_id_A = mol_id_A[0]
        self.mol_id_B = mol_id_B[0]
        self.group_id_A = mol_id_A[1]
        self.group_id_B = mol_id_B[1]
        self.true_dists_A = mol_id_A[2]
        self.true_dists_B = mol_id_B[2]
        self.real_site_A, self.real_site_group_id_A, \
            self.real_site_spe_A, self.real_site_mol_id_A = find_site(
                self.parent_A, self.mol_id_A, self.group_id_A, self.true_dists_A)
        self.real_site_B, self.real_site_group_id_B, \
            self.real_site_spe_B, self.real_site_mol_id_B = find_site(
                self.parent_B, self.mol_id_B, self.group_id_B, self.true_dists_B)
        self.mol_g_A = cal_g(struc_A, self.mol_id_A, self.group_id_A, self.true_dists_A)
        self.mol_g_B = cal_g(struc_B, self.mol_id_B, self.group_id_B, self.true_dists_B)
        self._rm_dummy = []
        self.mol_group_list_A = self.make_mol_group_list(self.real_site_mol_id_A, self.real_site_group_id_A)
        self.mol_group_list_B = self.make_mol_group_list(self.real_site_mol_id_B, self.real_site_group_id_B)
        count = 0
        # ---------- lattice crossover
        self._lattice_crossover()
        # ---------- generate children
        if len(self.parent_A) != len(self.parent_B):
            return None, None
        while True:
            count += 1
            # ------ coordinate crossover
            self._one_point_crossover_mol()
            self.child = Structure(lattice=self.lattice, species=self.species, coords=self.coords)
            self.decide_mol_from_dummy(self.child)

            # ------ check nat_diff
            self.check_mol_num()    # get self._nat_diff

            if any([abs(n) > rin.nat_diff_tole for n in self._nmol_diff]):
                if count > rin.maxcnt_ea:    # fail
                    self.child = None
                    return self.child, None
                continue    # slice again
            # ------ check mindist
            self.decide_mol_from_dummy(self.child)
            dist_list = self.check_distance_mol(self.child, rin.mindist_mol_ea, check_all=True)
            # ------ something smaller than mindist
            if dist_list:
                # -- remove atoms within mindist
                if any([n > 0 for n in self._nmol_diff]):
                    self._remove_within_mindist_mol()
                    if self.child is None:    # fail --> slice again
                        if count > rin.maxcnt_ea:
                            return None, None
                        continue
                else:    # nothing to remove, nat_diff = [0, 0]
                    if count > rin.maxcnt_ea:
                        return None, None
                    continue    # fail --> slice again
            # ------ recheck nmol_diff
            self.check_mol_num()
            # -- remove atoms near the border line
            if any([n > 0 for n in self._nmol_diff]):
                self._remove_border_line_mol()
            # -- add atoms near border line
            if any([n < 0 for n in self._nmol_diff]):
                self._add_border_line_mol()
            # -- success --> break while loop
            if self.child is not None:
                break
            # -- fail --> slice again
            else:
                if count > rin.maxcnt_ea:
                    return None, None
                continue
        # ---------- final check for nat
        self.decide_mol_from_dummy(self.child)
        self.child, self.child_mol_id, self.child_group_id = self.replace_mol(self.child)
        self._get_nat_diff()
        if not all([n == 0 for n in self._nat_diff]):
            return None, None    # failure
        # ---------- sort by atype
        self.child, self.child_mol_id, self.child_group_id = sort_by_atype_mol(self.child,
                                                                               rin.atype,
                                                                               self.child_mol_id,
                                                                               self.child_group_id)
        # ---------- return
        return self.child, [self.child_mol_id, self.child_group_id, self.true_dists_A]

    def _one_point_crossover_mol(self):
        # ---------- slice point
        while True:
            self._slice_point = np.random.normal(loc=0.5, scale=0.1)
            if 0.3 <= self._slice_point <= 0.7:
                break
        self._axis = np.random.choice([0, 1, 2])
        # ---------- crossover
        species_A = []
        species_B = []
        coords_A = []
        coords_B = []
        group_id_A = []
        group_id_B = []
        mol_id_A = []
        mol_id_B = []

        for i in range(len(self.mol_g_A)):
            if self.mol_g_A[i][self._axis] <= self._slice_point:
                species_A.append(DummySpecie("Xa{}".format(i)))
                coords_A.append(self.mol_g_A[i])
                group_id_A.append(i)
                for j, gid in enumerate(self.real_site_group_id_A):
                    if gid == i:
                        mol_id_A.append(self.real_site_mol_id_A[j])
                        break
            else:
                species_B.append(DummySpecie("Xa{}".format(i)))
                coords_B.append(self.mol_g_A[i])
                group_id_B.append(i)
                for j, gid in enumerate(self.real_site_group_id_A):
                    if gid == i:
                        mol_id_B.append(self.real_site_mol_id_A[j])
                        break
            if self.mol_g_B[i][self._axis] >= self._slice_point:
                species_A.append(DummySpecie("Xb{}".format(i)))
                coords_A.append(self.mol_g_B[i])
                group_id_A.append(i+len(self.group_id_A))
                for j, gid in enumerate(self.real_site_group_id_B):
                    if gid == i:
                        mol_id_A.append(self.real_site_mol_id_B[j])
                        break
            else:
                species_B.append(DummySpecie("Xb{}".format(i)))
                coords_B.append(self.mol_g_B[i])
                group_id_B.append(i+len(self.group_id_A))
                for j, gid in enumerate(self.real_site_group_id_B):
                    if gid == i:
                        mol_id_B.append(self.real_site_mol_id_B[j])
                        break
        # ---------- adopt a structure with more atoms
        if len(species_A) > len(species_B):
            species = species_A
            coords = coords_A
            group_id = group_id_A
            mol_id = mol_id_A

        elif len(species_A) < len(species_B):
            species = species_B
            coords = coords_B
            group_id = group_id_B
            mol_id = mol_id_B
        else:
            if np.random.choice([0, 1]):
                species = species_A
                coords = coords_A
                group_id = group_id_A
                mol_id = mol_id_A
            else:
                species = species_B
                coords = coords_B
                group_id = group_id_B
                mol_id = mol_id_B
        # ---------- set instance variables
        self.species, self.coords, self.child_mol_id, self.child_group_id = species, coords, mol_id, group_id

    def check_mol_num(self):
        self._nmol_diff = []
        for j, mol_number in enumerate(list(set(self.child_mol_id))):
            tmp_num = []
            for i, mol_id in enumerate(self.child_group_id):
                if mol_number == self.child_mol_id[i]:
                    tmp_num.append(mol_id)
            self._nmol_diff.append(len(list(set(tmp_num))) - rin.nmol[j])

    def _remove_border_line_mol(self):
        # ---------- rank atoms from border line
        coords_axis = self.child.frac_coords[:, self._axis]
        # ------ one point crossover: boundary --> 0.0, slice_point, 1.0
        near_sp = (self._slice_point/2.0 < coords_axis) & \
            (coords_axis < (self._slice_point + 1.0)/2.0)
        near_one = (self._slice_point + 1.0)/2.0 <= coords_axis
        # -- distance from nearest boundary
        coords_diff = np.where(near_sp,
                               abs(coords_axis - self._slice_point),
                               coords_axis)
        coords_diff = np.where(near_one, 1.0 - coords_diff, coords_diff)
        dummy_atom_border_indx = np.argsort(coords_diff)
        # ---------- remove list
        rm_list = []
        for itype, nrm in enumerate(self._nmol_diff):
            rm_list.append([])
            if nrm > 0:
                for dab_indx in dummy_atom_border_indx:
                    if self.child_mol_id[dab_indx] == itype:
                        rm_list[itype].append(dab_indx)
                    if len(rm_list[itype]) == nrm:
                        break
        for i, each_type in enumerate(rm_list):
            if each_type:
                self._nmol_diff[i] -= len(each_type)
        rm_list_one_division = sum(rm_list, [])
        # -- remove
        self.child.remove_sites(sorted(rm_list_one_division, reverse=True))
        for rl in sorted(rm_list_one_division, reverse=True):
            self.child_id.pop(rl)
            self.child_mol_id.pop(rl)
            self.child_group_id.pop(rl)
            self.child_parent.pop(rl)

    def _add_border_line_mol(self):
        for i in range(len(self._nmol_diff)):
            # ---------- counter
            cnt = 0
            # ---------- add atoms
            while self._nmol_diff[i] < 0:
                cnt += 1
                coords = np.random.rand(3)
                self._mean_choice()
                coords[self._axis] = np.random.normal(loc=self._mean,
                                                      scale=0.08)
                self.child.append(species=DummySpecie("Xx{}".format(i)), coords=coords)
                self._rm_dummy.append(DummySpecie("Xx{}".format(i)))
                self.decide_mol_from_dummy(self.child)
                success, mindist_ij, dist = self.check_distance_mol(self.child, rin.mindist_mol_ea)
                if success:
                    cnt = 0    # reset
                    self._nmol_diff[i] += 1
                else:
                    type0 = rin.atype[mindist_ij[0]]
                    type1 = rin.atype[mindist_ij[1]]
                    logger.warning(f'mindist in _add_border_line: {type0} - {type1}, {dist}. retry.')
                    self.child.pop()    # cancel
                # ------ fail
                if cnt == rin.maxcnt_ea:
                    self.child = None
                    return

    def replace_mol(self, struc):
        coords_AB = []
        species_AB = []
        mol_id_AB = []
        group_id_AB = []
        cnt = 0
        child_id = self.child_id
        child_parent = self.child_parent
        for i, dummy in enumerate(struc):
            if child_parent[i] == 0:
                for j, miA in enumerate(self.real_site_group_id_A):
                    if child_id[i] == miA:
                        tmp_coord = struc.cart_coords[i] + (self.parent_A.lattice.get_cartesian_coords(
                                                                self.real_site_A[j])
                                                            - self.parent_A.lattice.get_cartesian_coords(
                                                                self.mol_g_A[miA]))
                        coords_AB.append(struc.lattice.get_fractional_coords(tmp_coord))
                        species_AB.append(self.real_site_spe_A[j])
                        group_id_AB.append(cnt)
                        mol_id_AB.append(self.real_site_mol_id_A[j])
            if child_parent[i] == 1:
                for j, miB in enumerate(self.real_site_group_id_B):
                    if child_id[i] == miB:
                        tmp_coord = struc.cart_coords[i] + (self.parent_B.lattice.get_cartesian_coords(
                                                                self.real_site_B[j])
                                                            - self.parent_B.lattice.get_cartesian_coords(
                                                                self.mol_g_B[miB]))
                        coords_AB.append(struc.lattice.get_fractional_coords(tmp_coord))
                        species_AB.append(self.real_site_spe_B[j])
                        group_id_AB.append(cnt)
                        mol_id_AB.append(self.real_site_mol_id_B[j])
            cnt += 1
        for c_AB, s_AB in zip(coords_AB, species_AB):
            struc.append(s_AB, c_AB)
        struc.remove_species(self.species)
        struc.remove_species(self._rm_dummy)
        return struc, mol_id_AB, group_id_AB

    def decide_mol_from_dummy(self, struc):
        self.child_id = []
        self.child_parent = []
        for i, struc_spe in enumerate(struc.species):
            if struc_spe.symbol[1] == 'a':
                self.child_id.append(int(struc_spe.symbol[2:]))
                self.child_parent.append(0)
            elif struc_spe.symbol[1] == 'b':
                self.child_id.append(int(struc_spe.symbol[2:]))
                self.child_parent.append(1)
            elif struc_spe.symbol[1] == 'x':
                choice_parent = random.choice([0, 1])
                self.child_parent.append(choice_parent)
                if choice_parent:
                    choice_id = random.choice(self.mol_group_list_B[int(struc_spe.symbol[2:])])
                    self.child_id.append(choice_id)
                else:
                    choice_id = random.choice(self.mol_group_list_A[int(struc_spe.symbol[2:])])
                    self.child_id.append(choice_id)

    def check_distance_mol(self, struc, mindist_mol, check_all=False):
        # ---------- initialize
        if check_all:
            dist_list = []    # [(i, j, dist), (i, j, dist), ...]

        # ---------- in case there is only one atom
        if struc.num_sites == 1:
            dist = min(struc.lattice.abc)
            if dist < mindist_mol[0][0]:
                if check_all:
                    dist_list.append((0, 0, dist))
                    return dist_list
                return False, (0, 0), dist
            return True, None, None
        # ---------- normal case
        for i in range(struc.num_sites):
            for j in range(i):
                dist = struc.get_distance(j, i)
                if self.child_parent[i] == 0:
                    i_type = self.real_site_mol_id_A[i]
                else:
                    i_type = self.real_site_mol_id_B[i]
                if self.child_parent[j] == 0:
                    j_type = self.real_site_mol_id_A[j]
                else:
                    j_type = self.real_site_mol_id_B[j]
                if dist < mindist_mol[i_type][j_type]:
                    if check_all:
                        dist_list.append((j, i, dist))
                    else:
                        return False, (i_type, j_type), dist
        # ---------- return
        if check_all:
            if dist_list:
                dist_list.sort()    # sort
                return dist_list
            else:
                return dist_list    # dist_list is vacant list
        return True, None, None

    def make_mol_group_list(self, mol_id, group_id):
        cnt_id = 0
        mol_group_list = {}
        while True:
            if cnt_id not in mol_id:
                return mol_group_list
            tmp_group_id = []
            for i, mid in enumerate(mol_id):
                if cnt_id == mid:
                    tmp_group_id.append(group_id[i])
            mol_group_list.update({cnt_id: list(set(tmp_group_id))})
            cnt_id += 1

    def _remove_within_mindist_mol(self):
        '''
        if success: self.child <-- child structure data
        if fail: self.child <-- None
        '''
        for itype in range(len(rin.nmol)):
            while self._nmol_diff[itype] > 0:
                # ---------- check dist
                dist_list = self.check_distance_mol(self.child, rin.mindist_mol_ea, check_all=True)
                if not dist_list:    # nothing within mindist
                    return
                # ---------- appearance frequency
                ij_within_dist = [isite[0] for isite in dist_list] + [
                    jsite[1] for jsite in dist_list]
                site_counter = Counter(ij_within_dist)
                # ---------- get index for removing
                rm_index = []
                # ---- site[0]: index, site[1]: count
                for site in site_counter.most_common():
                    if self.child_parent[site[0]] == 0:
                        if self.real_site_mol_id_A[site[0]] == itype:
                            rm_index.append(site[0])
                    else:
                        if self.real_site_mol_id_B[site[0]] == itype:
                            rm_index.append(site[0])
                        break    # break for loop
                # ---------- remove atom
                if len(rm_index) == 0:
                    self.child = None
                    return
                else:
                    self.child.remove_sites(rm_index)
                    self._nmol_diff[itype] -= 1
                    sorted_rm_index = sorted(rm_index, reverse=True)
                    for indx in sorted_rm_index:
                        self.child_id.pop(indx)
                        self.child_parent.pop(indx)
                        self.child_group_id.pop(indx)
                        self.child_mol_id.pop(indx)
        # ---------- final check
        dist_list = self.check_distance_mol(self.child, rin.mindist_mol_ea, check_all=True)
        if dist_list:    # still something within mindist
            self.child = None
