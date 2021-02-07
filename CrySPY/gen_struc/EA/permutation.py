'''
Permutaion class
'''

import sys

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher

from ..struc_util import sort_by_atype, check_distance


class Permutation:
    '''
    permutation

    # ---------- args
    atype (list): atom type, e.g. ['Si', 'O'] for Si4O8

    mindist (2d list): constraint on minimum interatomic distance,
                       mindist must be a symmetric matrix
        e.g. [[1.8, 1.2], [1.2, 1.5]
            Si - Si: 1.8 angstrom
            Si - O: 1.2
            O - O: 1.5

    ntimes (int): ntimes permutaion

    maxcnt_ea (int): maximum number of trial in crossover
    # ---------- instance methods
    self.gen_child(struc)
        if success, return self.child
        if fail, return None
    '''

    def __init__(self, atype, mindist, ntimes=1, maxcnt_ea=100):
        # ---------- check args
        # ------ atype, nat, mindist
        for x in [atype, mindist]:
            if not isinstance(x, list):
                raise TypeError('{} must be list'.format(x))
        if not len(atype) == len(mindist):
            raise ValueError('not len(atype) == len(mindist)')
        # -- check symmetric
        for i in range(len(mindist)):
            for j in range(i):
                if not mindist[i][j] == mindist[j][i]:
                    raise ValueError('mindist is not symmetric. '
                                     '({}, {}): {}, ({}, {}): {}'.format(
                                         i, j, mindist[i][j],
                                         j, i, mindist[j][i]))
        self.atype = atype
        self.mindist = mindist
        # ------ ntimes, maxcnt_ea
        for x in [ntimes, maxcnt_ea]:
            if isinstance(x, int):
                if x <= 0:
                    raise ValueError('{} must be positive int'.format(x))
            else:
                raise TypeError('{} must be int'.format(x))
        self.ntimes = ntimes
        self.maxcnt_ea = maxcnt_ea

    def gen_child(self, struc):
        '''
        generate child struture

        # ---------- return
        (if success) self.child:
        (if fail) None:
        '''
        # ---------- keep original structure
        self.child = struc.copy()
        # ---------- instantiate StructureMatcher
        smatcher = StructureMatcher()    # instantiate
        # ---------- ntimes permutation
        cnt = 0
        while True:
            n = self.ntimes
            while n > 0:
                # ------ prepare index for each atom type
                indx_each_type = []
                for a in self.atype:
                    indx_each_type.append(
                        [i for i, site in enumerate(self.child)
                         if site.species_string == a])
                # ------ choose two atom type
                type_choice = np.random.choice(len(self.atype), 2,
                                               replace=False)
                # ------ choose index
                indx_choice = []
                for tc in type_choice:
                    indx_choice.append(np.random.choice(indx_each_type[tc]))
                # ------ replace each other
                self.child.replace(indx_choice[0],
                                   species=self.atype[type_choice[1]])
                self.child.replace(indx_choice[1],
                                   species=self.atype[type_choice[0]])
                # ------ compare to original one
                if smatcher.fit(self.child, struc):
                    n = self.ntimes    # back to the start
                    continue
                else:
                    n -= 1
            # ------ check distance
            success, mindist_ij, dist = check_distance(self.child,
                                                       self.atype,
                                                       self.mindist)
            if success:
                self.child = sort_by_atype(self.child, self.atype)
                return self.child
            else:
                sys.stderr.write('mindist in permutation: {} - {}, {}. retry.\n'.format(
                    self.atype[mindist_ij[0]],
                    self.atype[mindist_ij[1]],
                    dist))
                cnt += 1
                if cnt >= self.maxcnt_ea:
                    print('Permutatin: could not satisfy min_dist' +
                          ' in {} times'.format(self.maxcnt_ea))
                    print('Change parent')
                    self.child = None
                    return None    # change parent
