#!/usr/bin/env python
# -*- coding: utf-8 -*-
# --------------------------------------------------------------------------
#
# This code partly includes find_wy (https://github.com/nim-hrkn/find_wy)
# which is distributed under the Apache License, Version 2.0.
#
# --------------------------------------------------------------------------

from __future__ import print_function

import json

import numpy as np
from pymatgen import Structure

from ..dist import check_min_dist


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
    n_uniq, wydata_eq_atom = get_wydata_eq_atom(wydata)
    eq_atomnames = {}
    eq_positions = {}
    for key, value in sorted(n_uniq.items(), key=lambda x: x[1]):    # equivalent atom loop
        # ------ distribute eq atoms. first, special (num_uniqvar = 0), then, others
        cnt = 0
        while True:
            eq_atomnames[key], eq_positions[key] = gen_eq_atoms(wydata_eq_atom[key])
            # -- sort in original order
            atomnames = []
            positions = []
            for key_a, value_a in sorted(eq_atomnames.items()):
                atomnames += eq_atomnames[key_a]
                positions += eq_positions[key_a]
            # -- Cartesian coordinate; use clat (not plat)
            cart = []
            for p in positions:
                v = np.zeros(3)
                for i in range(3):
                    a = np.array(clat[i])
                    v += p[i] * a
                cart.append(v)
            # -- check minimum distance
            spgstruc = Structure(plat, atomnames, cart, coords_are_cartesian=True)
            if check_min_dist(spgstruc, cumul_nat, mindist) is False:
                cnt = maxcnt + 1 if value == 0 else cnt + 1    # num_uniqvar = 0 --> value == 0
                if maxcnt < cnt:
                    return False, spgstruc    # spgstruc is dummy
            else:
                break    # break while loop --> next eq atoms
    return True, spgstruc


def get_wydata_eq_atom(wydata):
    i = 0    # count eq_atom, not atom
    n_uniq = {}    # num_uniqvar each eq_atom
    wydata_eq_atom = {}    # wydata each eq_atom
    for specie in wydata['atoms']:
        for wydata2 in specie:    # equivalent atom loop
            n_uniq[i] = wydata2[0]['num_uniqvar']
            wydata_eq_atom[i] = wydata2
            i += 1
    return n_uniq, wydata_eq_atom


def gen_eq_atoms(wydata2):
    eq_atomnames = []
    eq_positions = []
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
        eq_positions.append(pos + each['add'])
        eq_atomnames.append(each['name'])

    return eq_atomnames, eq_positions
