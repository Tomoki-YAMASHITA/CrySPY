
from logging import getLogger

import numpy as np

from ...IO import read_input as rin
from ...util.struc_util import check_distance, sort_by_atype
#from .adj_comp import operation_atoms, convex_hull_check


logger = getLogger('cryspy')


class Substitution:

    def __init__(self, mindist, target='random'):
        self.mindist = mindist
        self.target = target

    def gen_child(self, struc, atype_avail_add, atype_avail_elim):
        cnt = 0
        while True:
            # ---------- keep original structure
            self.child = struc.copy()
            # ---------- substitution
            if self.target == 'random':
                # ------ random choice for atom type
                at = np.random.choice(atype_avail_elim)
                if at in atype_avail_add:
                    atype_avail_add.remove(at)
                if len(atype_avail_add) == 0:
                    logger.warning('Substitution: no atype_avail_add, retry.')
                    cnt += 1
                    if cnt >= rin.maxcnt_ea:
                        logger.warning('Substitution: cnt >= maxcnt_ea.')
                        logger.warning('Change parent')
                        self.child = None
                        return None
                    continue
                else:
                    # ------ elimination
                    aindx = [i for i, site in
                             enumerate(self.child) if site.species_string == at]
                    elim_indx = np.random.choice(aindx, 1)  # ", 1" to get array
                    coords = self.child[elim_indx[0]].frac_coords    # coords of eliminated atom
                    self.child.remove_sites(elim_indx)
                    # ------ addition
                    at = np.random.choice(atype_avail_add)
                    self.child.append(species=at, coords=coords)
            # ---------- check mindist
            success, mindist_ij, dist = check_distance(self.child,
                                                       rin.atype,
                                                       self.mindist)
            if success:
                self.child = sort_by_atype(self.child, rin.atype)
                return self.child
            else:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist in addition: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                if cnt >= rin.maxcnt_ea:
                    logger.warning('Substitution: cnt >= maxcnt_ea.')
                    logger.warning('Change parent')
                    self.child = None
                    return None
