'''
Rotation class
'''

from logging import getLogger
import sys

import numpy as np
from pymatgen.core import Structure

from ...IO import read_input as rin
from ...util.struc_util import check_distance, rot_mat, find_site, cal_g, sort_by_atype_mol


logger = getLogger('cryspy')

class Rotation:
    '''
    rotation
    '''

    def __init__(self, mindist):
        self.mindist = mindist

    def gen_child(self, struc, mol_id):
        # ---------- init
        cnt = 0
        self.child = struc
        real_coords, real_group_id, real_species, _ = find_site(self.child, mol_id[0], mol_id[1], mol_id[2])
        mol_gs = self.child.lattice.get_cartesian_coords(cal_g(self.child, mol_id[0], mol_id[1], mol_id[2]))
        rotated_group_id = []
        # ---------- generate rotated structure
        while True:
            cnt += 1
            rot_coords = []
            rot_species = []
            # ------ rotation molecules
            for i in range(sum(rin.nmol)):
                # -- calculate rotation matrix
                mol_angle = np.deg2rad(rin.rot_max_angle) * np.random.rand(3)
                R = rot_mat(mol_angle)
                # -- rotate
                for j in range(len(self.child)):
                    if real_group_id[j] == i:
                        rot_coords.append(np.matmul(self.child.lattice.get_cartesian_coords(real_coords[j]) -
                                                    mol_gs[i], R) + mol_gs[i])
                        rot_species.append(real_species[j])
                        rotated_group_id.append(i)
            # ------ child
            self.child = Structure(lattice=self.child.lattice, species=rot_species, coords=rot_coords,
                                   coords_are_cartesian=True)

            # ------ check distance
            success, mindist_ij, dist = check_distance(self.child,
                                                       rin.atype,
                                                       self.mindist)
            if success:
                self.child, self.mol_id, self.group_id = sort_by_atype_mol(self.child, rin.atype,
                                                                           mol_id[0], rotated_group_id)
                return self.child, [self.mol_id, self.group_id, mol_id[2]]
            else:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist in rotation: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                if cnt >= rin.maxcnt_ea:
                    logger.warning('Rotation: could not satisfy min_dist' +
                          f' in {rin.maxcnt_ea} times')
                    logger.warning('Change parent')
                    self.child = None
                    return None, None    # change parent
