from logging import getLogger

import numpy as np

from ...IO import read_input as rin
from ...util.struc_util import check_distance, sort_by_atype
#from .adj_comp import operation_atoms, convex_hull_check


logger = getLogger('cryspy')


class Addition:

    def __init__(self, mindist, target='random'):
        self.mindist = mindist
        self.target = target

    def gen_child(self, struc, atype_avail):
        cnt = 0
        while True:
            # ---------- keep original structure
            self.child = struc.copy()
            # ---------- addition
            if self.target == 'random':
                # ------ random choice for atom type
                at = np.random.choice(atype_avail)
                # ------ add atom
                coords = np.random.rand(3)
                self.child.append(species=at, coords=coords)
            # ---------- not implemented yet
            # elif target in ['depop', 'overpop']:
            #     section = convex_hull_check()
            #     self.child = operation_atoms('addition', self.child, section)
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
                logger.warning(f'mindist in addition: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                if cnt >= rin.maxcnt_ea:
                    logger.warning('Addition: could not satisfy min_dist' +
                          f' in {rin.maxcnt_ea} times')
                    logger.warning('Change parent')
                    self.child = None
                    return None
