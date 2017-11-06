#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function


def check_min_dist(structure, cumul_nat, mindist):
    if structure.num_sites == 1:
        return True
    for i in xrange(structure.num_sites):
        for j in xrange(structure.num_sites):
            if i < j:
                dist = structure.get_distance(i, j)
                type_i = get_atype_num(cumul_nat, i)
                type_j = get_atype_num(cumul_nat, j)
                if dist < mindist[type_i][type_j]:
                    return False
    return True


def get_atype_num(cumul_nat, site_num):
    '''
    e.g. SrTiO3
    atype = ['Sr', 'Ti', 'O']
    nat = [1, 1, 3]

    atom 0: Sr1 --> atype_num = 0
    atom 1: Ti1 --> atype_num = 1
    atom 2: O1 --> atype_num = 2
    atom 3: O2 --> atype_num = 2
    atom 4: O3 --> atype_num = 2
    '''
    for atype_num, cumul_i in enumerate(cumul_nat):
        if site_num < cumul_i:
            return atype_num
