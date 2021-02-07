'''
Strain class
'''
import sys

import numpy as np
from pymatgen import Structure

from ..struc_util import sort_by_atype, check_distance


class Strain:
    '''
    strain

    # ---------- args
    sigma (float): standard deviation for strain

    # ---------- instance methods
    gen_child(self, struc)
    '''

    def __init__(self, atype, mindist, sigma=0.5, maxcnt_ea=100):
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
        # ------ sigma
        if isinstance(sigma, float):
            if sigma <= 0:
                raise ValueError('simga must be positive float')
        else:
            raise TypeError('sigma must be float')
        self.sigma = sigma
        # ------ maxcnt_ea
        if isinstance(maxcnt_ea, int):
            if maxcnt_ea <= 0:
                raise ValueError('maxcnt_ea must be positive int')
        else:
            raise TypeError('maxcnt_ea must be int')
        self.maxcnt_ea = maxcnt_ea

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
                                loc=0.0, scale=self.sigma)
                        else:
                            strain_matrix[i][j] = np.random.normal(
                                loc=0.0, scale=self.sigma)/2.0
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
                    self.child = None
                    return None    # change parent
