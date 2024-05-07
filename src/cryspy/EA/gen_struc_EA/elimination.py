from logging import getLogger

import numpy as np
from pymatgen.core import Structure

from ...IO import read_input as rin
from ...util.struc_util import sort_by_atype
#from .adj_comp import operation_atoms, convex_hull_check


logger = getLogger('cryspy')

class Elimination:
    
    def __init__(self, mindist, target='random'):
        self.mindist = mindist
        self.target = target

    def gen_child(self, struc, atype_avail):
        # ---------- keep original structure
        self.child = struc.copy()
        # ---------- elimination
        if self.target == 'random':
            # ------ random choice for atom type
            at = np.random.choice(atype_avail)
            # ------ random choice for atom index
            aindx = [i for i, site in enumerate(self.child) if site.species_string == at]
            elim_indx = np.random.choice(aindx, 1)  # ", 1" to get array
            # ------ remove atom
            self.child.remove_sites(elim_indx)
        # ---------- not implemented yet
        # elif tgt in ['depop', 'overpop']:
        #     section = convex_hull_check()
        #     self.child = operation_atoms('elimination', self.child, section)
        # ---------- return
        #            no need to check distance in elimination
        return self.child
