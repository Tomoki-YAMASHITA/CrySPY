#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from pymatgen.io.cif import CifWriter

from ..IO import read_input as rin


def out_opt_struc(opt_struc, current_id):
    # ---------- cut unnecessary parts in poscar
    pos = opt_struc.to(fmt='poscar')
    pos = pos.split('\n')
    blank_indx = pos.index('')
    pos = pos[:blank_indx]
    pos[0] = 'ID_{}'.format(current_id)    # replace with ID
    lines = [line+'\n' for line in pos]

    # ---------- add in ./data/opt_POSCARS
    with open('./data/opt_POSCARS', 'a') as fopt:
        for line in lines:
            fopt.write(line)


def out_opt_cif(opt_struc, current_id, work_path):
    # ---------- opt_CIFS
    cif = CifWriter(opt_struc, symprec=rin.symtoleR)
    cif.write_file(work_path+'tmp.cif')

    # ---------- correct title (need to delete '_chemical_formula_sum')
    with open(work_path+'tmp.cif', 'r') as fcif:
        ciflines = fcif.readlines()
    ciflines[1] = 'data_ID_{}\n'.format(current_id)
    if ciflines[11][:21] == '_chemical_formula_sum':
        ciflines.pop(11)
    else:
        raise ValueError('ciflines[11] is not _chemical_formula_sum, have to fix bag')

    # ---------- cif --> opt_cifs
    with open('./data/opt_CIFS.cif', 'a') as foptcif:
        for line in ciflines:
            foptcif.write(line)

    # ---------- clean tmp.cif
    os.remove(work_path+'tmp.cif')
