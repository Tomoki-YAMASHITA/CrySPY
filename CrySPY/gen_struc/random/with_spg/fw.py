#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------
#
# This code partly includes find_wy (https://github.com/nim-hrkn/find_wy)
# which is distributed under the Apache License, Version 2.0.
#
# --------------------------------------------------------------------------

from __future__ import print_function

import copy
import json

import numpy as np
from pymatgen import Structure


def fw_input(atype, nat, spg, a, b, c, cosa, cosb, cosg):
    with open('input', 'w') as f:
        f.write('nspecies {}\n'.format(len(atype)))
        f.write('species_name')
        for aa in atype:
            f.write('  {}'.format(aa))
        f.write('\n')
        f.write('species_num')
        for i in nat:
            f.write('  {}'.format(i))
        f.write('\n')
        f.write('spacegroup  {}\n'.format(spg))
        f.write('originchoice  1\n')
        f.write('\n')
        f.write('a  {}\n'.format(a))
        f.write('b  {}\n'.format(b))
        f.write('c  {}\n'.format(c))
        f.write('cosa  {}\n'.format(cosa))
        f.write('cosb  {}\n'.format(cosb))
        f.write('cosc  {}\n'.format(cosg))
        f.write('\n')
        # f.write('selectone true\n')
        f.write('randomseed auto\n')


def gen_wypos(cumul_nat, mindist, maxcnt):
    '''
    Success --> return True, structure data
    Failure --> return False, _
    '''

    # ---------- load POS_WY_SKEL_ALL.json
    with open('POS_WY_SKEL_ALL.json', 'r') as f:
        wydata = json.load(f)

    # ---------- generate structure
    plat = wydata['primitivevector']
    clat = wydata['conventionalvector']
    atomnames = []
    positions = []
    for specie in wydata['atoms']:
        for wydata2 in specie:    # equivalent atom loop
            cnt = 0
            while True:
                tmp_atomnames, tmp_positions = gen_eq_atoms(wydata2, atomnames, positions)

                # ------ Cartesian coordinate
                #      platではなくclatを使って変換しないと上手くいかない
                cart = []
                for p in tmp_positions:
                    v = np.zeros(3)
                    for i in range(3):
                        a = np.array(clat[i])
                        v += p[i] * a
                    cart.append(v)

                # ------ check minimum distance
                spgstruc = Structure(plat, tmp_atomnames, cart, coords_are_cartesian=True)
                if check_min_dist(spgstruc, cumul_nat, mindist) is False:
                    cnt += 1
                    if maxcnt < cnt:
                        return False, spgstruc    # spgstruc is dummy
                else:
                    atomnames, positions = tmp_atomnames, tmp_positions
                    break
    return True, spgstruc


def gen_eq_atoms(wydata2, atomnames, positions):
    tmp_atomnames = copy.deepcopy(atomnames)
    tmp_positions = copy.deepcopy(positions)
    rval = np.random.random_sample(3)
    for each in wydata2:
        pos = []
        for ch in each['xyzch']:
            if ch == '-2x':
                pos.append(-2.0 * rval[0])
            elif ch == '-x+y':
                pos.append(-rval[0] + rval[1])
            elif ch == '-z':
                pos.append(-rval[2])
            elif ch == '-y':
                pos.append(-rval[1])
            elif ch == '-x':
                pos.append(-rval[0])
            elif ch == '0':
                pos.append(0.0)
            elif ch == 'x':
                pos.append(rval[0])
            elif ch == 'y':
                pos.append(rval[1])
            elif ch == 'z':
                pos.append(rval[2])
            elif ch == 'x-y':
                pos.append(rval[0] - rval[1])
            elif ch == '2x':
                pos.append(2.0 * rval[0])
            else:
                raise ValueError('unknown ch in conversion in gen_wycoord')
        pos = np.array(pos)
        tmp_positions.append(pos + each['add'])
        tmp_atomnames.append(each['name'])

    return tmp_atomnames, tmp_positions


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


# ---------- obsolete
#                version 0.4.x or lower
# def check_min_dist(structure):
#     if structure.num_sites == 1:
#         return 100.0    # dummy

#     min_dist = structure.get_distance(0, 1)
#     for i in xrange(structure.num_sites):
#         for j in xrange(structure.num_sites):
#             if i < j:
#                 dist = structure.get_distance(i, j)
#                 if dist < min_dist:
#                     min_dist = dist
#     return min_dist


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
