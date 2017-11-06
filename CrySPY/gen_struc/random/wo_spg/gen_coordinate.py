#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import random

from pymatgen import Structure

from ..dist import check_min_dist


def rndgen_coord(natot, va, vb, vc, atomlist, cumul_nat, mindist, maxcnt):
    '''
    Success --> return structure data in pymatgen format
    Failure --> return None
    '''
    # ---------- generate internal coordinates
    cnt = 0
    incoord = []
    while len(incoord) < natot:
        tmp_coord = [random.random() for i in range(3)]
        incoord.append(tmp_coord)
        tmp_struc = Structure([va, vb, vc], atomlist[:len(incoord)], incoord)
        if check_min_dist(tmp_struc, cumul_nat, mindist) is False:
            incoord.pop(-1)    # cancel
            cnt += 1
            if maxcnt < cnt:
                return None
    return tmp_struc
