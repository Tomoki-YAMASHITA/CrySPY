#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from pymatgen import Structure

from ...struc_util import check_distance


def rndgen_coord(natot, atype, nat, va, vb, vc, mindist, maxcnt):
    '''
    Success --> return structure data in pymatgen format
    Failure --> return None
    '''
    # ---------- initialize
    cnt = 0
    incoord = []
    atomlist = get_atomlist(atype, nat)

    # ---------- generate internal coordinates
    while len(incoord) < natot:
        tmp_coord = np.random.rand(3)
        incoord.append(tmp_coord)
        tmp_struc = Structure([va, vb, vc], atomlist[:len(incoord)], incoord)
        if not check_distance(tmp_struc, atype, mindist):
            incoord.pop(-1)    # cancel
            cnt += 1
            if maxcnt < cnt:
                return None
    return tmp_struc


def get_atomlist(atype, nat):
    '''
    e.g. Na2Cl2
        atomlist = ['Na', 'Na', 'Cl', 'Cl']
    '''
    atomlist = []
    for i in range(len(atype)):
        atomlist += [atype[i]]*nat[i]
    return atomlist
