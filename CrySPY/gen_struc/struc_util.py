#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np
from pymatgen import Structure
from pymatgen.io.cif import CifWriter


def out_poscar(struc, cID, fpath):
    # ---------- poscar format
    pos = struc.to(fmt='poscar')
    pos = pos.split('\n')
    blank_indx = pos.index('')    # cut unnecessary parts
    pos = pos[:blank_indx]
    pos[0] = 'ID_{}'.format(cID)    # replace with ID
    lines = [line+'\n' for line in pos]

    # ---------- append POSCAR
    with open(fpath, 'a') as f:
        for line in lines:
            f.write(line)


def out_cif(struc, cID, tmp_path, fpath, symprec=0.1):
    # ---------- opt_CIFS
    cif = CifWriter(struc, symprec=symprec)
    cif.write_file(tmp_path+'tmp.cif')

    # ---------- correct title for VESTA (need to delete '_chemical_formula_sum')
    with open(tmp_path+'tmp.cif', 'r') as fcif:
        ciflines = fcif.readlines()
    ciflines[1] = 'data_ID_{}\n'.format(cID)
    if ciflines[11][:21] == '_chemical_formula_sum':
        ciflines.pop(11)
    else:
        raise ValueError('ciflines[11] is not _chemical_formula_sum, have to fix bag')

    # ---------- cif --> opt_cifs
    with open(fpath, 'a') as foptcif:
        for line in ciflines:
            foptcif.write(line)

    # ---------- clean tmp.cif
    os.remove(tmp_path+'tmp.cif')


def frac_coord_zero_one(struc_in):
    '''
    fractional coordinates: 0.0 <= x,y,z < 1.0
    e.g. [0.0, -0.25, 0.7] --> [0.0, 0.75, 0.7]

    # ---------- args
    struc_in: structure data in pymatgen format
    '''
    struc = struc_in.copy()
    for i in range(struc.num_sites):
        struc[i] = struc[i].to_unit_cell
    return struc


def origin_shift(struc_in):
    '''
    Randomly shift the origin of struc_in

    # ---------- args
    struc_in: structure data in pymatgen format

    # ---------- return
    origin shifted structure (not change original struc_in)
    '''
    struc = struc_in.copy()
    coords_trans = struc.frac_coords + np.random.rand(3)
    struc_shift = Structure(struc.lattice, struc.species, coords_trans)
    struc_shift = frac_coord_zero_one(struc_shift)
    return struc_shift


def sort_by_atype(struc, atype):
    '''
    return a structre sorted by atype order as a new structure
    '''
    return struc.get_sorted_structure(key=lambda x: atype.index(x.species_string))


def check_distance(struc, atype, mindist, check_all=False):
    '''
    # ---------- args
    struc: structure data in pymatgen format
    atype (list): e.g. ['Li', 'Co, 'O']
    mindist (2d list) : e.g. [[2.0, 2.0, 1.2], [2.0, 2.0, 1.2], [1.2, 1.2, 1.5]]
    check_all (bool) : if True, check all atom pairs, and return True or False, dist_list
                       if False, stop when smaller distance (dist < mindist) is found
    # ---------- return
    (check_all=False) True: nothing smaller than mindist
    (check_all=False) False: something smaller than mindst
    (check_all=True) dist_list: if dist_list is vacant, nothing smaller than mindist
    '''

    # ---------- initialize
    if check_all:
        dist_list = []    # [(i, j, dist), (i, j, dist), ...]

    # ---------- in case there is only one atom
    if struc.num_sites == 1:
        dist = min(struc.lattice.abc)
        if dist < mindist[0][0]:
            if check_all:
                dist_list.append((0, 0, dist))
                return dist_list
            return False
        return True

    # ---------- normal case
    for i in range(struc.num_sites):
        for j in range(i):
            dist = struc.get_distance(j, i)
            i_type = atype.index(struc[i].species_string)
            j_type = atype.index(struc[j].species_string)
            if dist < mindist[i_type][j_type]:
                if check_all:
                    dist_list.append((j, i, dist))
                else:
                    return False

    # ---------- return
    if check_all:
        if dist_list:
            dist_list.sort()    # sort
            return dist_list
        else:
            return dist_list    # dist_list is vacant list
    return True
