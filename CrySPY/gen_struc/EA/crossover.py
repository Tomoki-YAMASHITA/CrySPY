'''
Crossover class
'''

from collections import Counter
import sys

import numpy as np
from pymatgen import Structure, Lattice

from ..struc_util import origin_shift, sort_by_atype, check_distance


class Crossover:
    '''
    crossover

    # ---------- args
    atype (list): atom type, e.g. ['Si', 'O'] for Si4O8

    nat (list): number of atom, e.g. [4, 8] for Si4O8

    mindist (2d list): constraint on minimum interatomic distance,
                       mindist must be a symmetric matrix
        e.g. [[1.8, 1.2], [1.2, 1.5]
            Si - Si: 1.8 angstrom
            Si - O: 1.2
            O - O: 1.5

    crs_lat ('equal' or 'random') how to mix lattice vectors

    nat_diff_tole (int): tolerance for difference in number of atoms
                         in crossover

    maxcnt_ea (int): maximum number of trial in crossover

    # ---------- instance methods
    self.gen_child(struc_A, struc_B)
        if success, return self.child
        if fail, return None
    '''

    def __init__(self, atype, nat, mindist, crs_lat='equal',
                 nat_diff_tole=4, maxcnt_ea=100):
        # ---------- check args
        # ------ atype, nat, mindist
        for x in [atype, nat, mindist]:
            if type(x) is not list:
                raise ValueError('atype, nat, and mindist must be list')
        if not len(atype) == len(nat) == len(mindist):
            raise ValueError('not len(atype) == len(nat) == len(mindist)')
        # -- check symmetric
        for i in range(len(mindist)):
            for j in range(i):
                if not mindist[i][j] == mindist[j][i]:
                    raise ValueError('mindist is not symmetric. ' +
                                     '({}, {}): {}, ({}, {}): {}'.format(
                                         i, j, mindist[i][j],
                                         j, i, mindist[j][i]))
        self.atype = atype
        self.nat = nat
        self.mindist = mindist
        # ------ crs_lat
        if crs_lat == 'equal':
            self.w_lat = np.array([1.0, 1.0])
        elif crs_lat == 'random':
            self.w_lat = np.random.choice([0.0, 1.0], size=2, replace=False)
        else:
            raise ValueError('crs_lat must be equal or random')
        # ------ nat_diff_tole, maxcnt_ea
        for x in [nat_diff_tole, maxcnt_ea]:
            if type(x) is int and x > 0:
                pass
            else:
                raise ValueError('nat_diff_tole and maxcnt_ea'
                                 ' must be positive int')
        self.nat_diff_tole = nat_diff_tole
        self.maxcnt_ea = maxcnt_ea

    def gen_child(self, struc_A, struc_B):
        '''
        generate child struture

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
            self._check_nat()    # get self._nat_diff
            if any([abs(n) > self.nat_diff_tole for n in self._nat_diff]):
                if count > self.maxcnt_ea:    # fail
                    self.child = None
                    return self.child
                continue    # slice again
            # ------ check mindist
            dist_list = check_distance(self.child, self.atype,
                                       self.mindist, check_all=True)
            # ------ something smaller than mindist
            if dist_list:
                # -- remove atoms within mindist
                if any([n > 0 for n in self._nat_diff]):
                    self._remove_within_mindist()
                    if self.child is None:    # fail --> slice again
                        if count > self.maxcnt_ea:
                            return None
                        continue
                else:    # nothing to remove, nat_diff = [0, 0]
                    if count > self.maxcnt_ea:
                        return None
                    continue    # fail --> slice again
            # ------ recheck nat_diff
            self._check_nat()
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
                if count > self.maxcnt_ea:
                    return None
                continue
        # ---------- final check for nat
        self._check_nat()
        if not all([n == 0 for n in self._nat_diff]):
            return None    # failure
        # ---------- sort by atype
        self.child = sort_by_atype(self.child, self.atype)
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

    def _check_nat(self):
        self._nat_diff = []
        species_list = [a.species_string for a in self.child]
        for i in range(len(self.atype)):
            self._nat_diff.append(species_list.count(self.atype[i])
                                  - self.nat[i])

    def _remove_within_mindist(self):
        '''
        if success: self.child <-- child structure data
        if fail: self.child <-- None
        '''
        for itype in range(len(self.atype)):
            while self._nat_diff[itype] > 0:
                # ---------- check dist
                dist_list = check_distance(self.child, self.atype,
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
                    if self.child[site[0]].species_string == self.atype[itype]:
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
        dist_list = check_distance(self.child, self.atype,
                                   self.mindist, check_all=True)
        if dist_list:    # still something within mindist
            self.child = None

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
                    if self.child[ab_indx].species_string == self.atype[itype]:
                        rm_list[itype].append(ab_indx)
                    if len(rm_list[itype]) == nrm:
                        break
        # ---------- remove
        for each_type in rm_list:
            if each_type:
                self.child.remove_sites(each_type)

    def _add_border_line(self):
        for i in range(len(self.atype)):
            # ---------- counter
            cnt = 0
            # ---------- add atoms
            while self._nat_diff[i] < 0:
                cnt += 1
                coords = np.random.rand(3)
                self._mean_choice()
                coords[self._axis] = np.random.normal(loc=self._mean,
                                                      scale=0.08)
                self.child.append(species=self.atype[i], coords=coords)
                success, mindist_ij, dist = check_distance(self.child,
                                                           self.atype,
                                                           self.mindist)
                if success:
                    cnt = 0    # reset
                    self._nat_diff[i] += 1
                else:
                    sys.stderr.write('mindist in _add_border_line: {} - {}, {}. retry.\n'.format(
                        self.atype[mindist_ij[0]],
                        self.atype[mindist_ij[1]],
                        dist))
                    self.child.pop()    # cancel
                # ------ fail
                if cnt == self.maxcnt_ea:
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
