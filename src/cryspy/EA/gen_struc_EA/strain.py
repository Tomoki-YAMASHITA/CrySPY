'''
Strain class
'''
from logging import getLogger
import sys

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie

from ...IO import read_input as rin
from ...util.struc_util import sort_by_atype, sort_by_atype_mol, check_distance, cal_g, find_site


logger = getLogger('cryspy')

class Strain:
    '''
    strain

    # ---------- instance methods
    gen_child(self, struc)
    '''

    def __init__(self, mindist):
        self.mindist = mindist

    def gen_child(self, struc):
        '''
        generate child struture

        # ---------- return
        (if success) self.child:
        (if fail) None:
        '''
        # ---------- init
        cnt = 0
        lat_mat = struc.lattice.matrix.T    # lattice vector as matrix
        # ---------- generate strained structure
        while True:
            # ------ strain matrix
            strain_matrix = np.empty([3, 3])
            for i in range(3):
                for j in range(3):
                    if i <= j:
                        if i == j:
                            strain_matrix[i][j] = 1.0 + np.random.normal(
                                loc=0.0, scale=rin.sigma_st)
                        else:
                            strain_matrix[i][j] = np.random.normal(
                                loc=0.0, scale=rin.sigma_st)/2.0
                            strain_matrix[j][i] = strain_matrix[i][j]
            # ------ strained lattice
            strained_lattice = np.dot(strain_matrix, lat_mat).T
            # ------ child
            self.child = Structure(strained_lattice, struc.species,
                                   struc.frac_coords)
            # ------ scale lattice
            self.child.scale_lattice(struc.volume)
            # ------ check distance
            success, mindist_ij, dist = check_distance(self.child,
                                                       rin.atype,
                                                       self.mindist)
            if success:
                self.child = sort_by_atype(self.child, rin.atype)
                return self.child
            else:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist in strain: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                if cnt >= rin.maxcnt_ea:
                    self.child = None
                    return None    # change parent

    def gen_child_mol(self, struc, mol_id):
        '''
        generate child structures for mol

        # ---------- return
        (if success) self.child:
        (if fail) None:
        '''
        # ---------- init
        cnt = 0
        lat_mat = struc.lattice.matrix.T    # lattice vector as matrix
        # ---------- generate strained structure
        while True:
            # ------ strain matrix
            strain_matrix = np.empty([3, 3])
            for i in range(3):
                for j in range(3):
                    if i <= j:
                        if i == j:
                            strain_matrix[i][j] = 1.0 + np.random.normal(
                                loc=0.0, scale=rin.sigma_st)
                        else:
                            strain_matrix[i][j] = np.random.normal(
                                loc=0.0, scale=rin.sigma_st)/2.0
                            strain_matrix[j][i] = strain_matrix[i][j]
            # ------ strained lattice
            strained_lattice = np.dot(strain_matrix, lat_mat).T
            # -- generate dummy structure
            dum_coords = cal_g(struc, mol_id[0], mol_id[1], mol_id[2])
            dum_species = []
            for i in range(len(dum_coords)):
                dum_species.append(DummySpecie("X{}".format(i)))
            dum_struc = Structure(strained_lattice, dum_species, dum_coords)
            # ------ scale lattice
            dum_struc.scale_lattice(struc.volume)
            # --calcurate atom coords to replace molecure
            fix_frac_coords, fix_group_id, fix_species, fix_mol_id = find_site(struc,
                                                                               mol_id[0],
                                                                               mol_id[1],
                                                                               mol_id[2])
            rm_spe = []
            strained_mol_id = []
            strained_group_id = []
            append_spe_co = []
            append_mol_gr = []
            for i, dum in enumerate(dum_struc):
                for j, mid in enumerate(fix_group_id):
                    if mid == i:
                        coord = (struc.lattice.get_cartesian_coords(fix_frac_coords[j])
                                 - struc.lattice.get_cartesian_coords(dum_coords[i])
                                 + dum_struc.cart_coords[i])
                        append_spe_co.append([fix_species[j], coord])
                        append_mol_gr.append([fix_mol_id[j], fix_group_id[j]])
                rm_spe.append(dum.specie)
            # replace mol
            for spe_co, mol_gr in zip(append_spe_co, append_mol_gr):
                dum_struc.append(spe_co[0], spe_co[1], coords_are_cartesian=True)
                strained_mol_id.append(mol_gr[0])
                strained_group_id.append(mol_gr[1])
            # remove dummy
            dum_struc.remove_species(rm_spe)
            # ------ child
            self.child = dum_struc
            # ------ check distance
            success, mindist_ij, dist = check_distance(self.child,
                                                       rin.atype,
                                                       self.mindist)
            if success:
                self.child, self.mol_id, self.group_id = sort_by_atype_mol(self.child, rin.atype,
                                                                           strained_mol_id, strained_group_id)
                return self.child, [self.mol_id, self.group_id, mol_id[2]]
            else:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist in permutation: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                if cnt >= rin.maxcnt_ea:
                    self.child = None
                    return None    # change parent
