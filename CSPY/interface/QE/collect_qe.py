#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import numpy as np
from pymatgen.io.cif import CifWriter

from . import structure as qe_structure
from ...IO import read_input as rin


def collect_qe(current_id, work_path):
    #---------- check optimization in previous stage
    try:
        with open(work_path+'prev_'+rin.outfile, 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'End final coordinates' in line:
                check_opt = 'done'
    except:
        check_opt = 'no_file'

    #---------- obtain energy and magmom
    try:
        with open(work_path+rin.outfile, 'r') as fpout:
            lines = fpout.readlines()
        energy = np.nan
        for line in lines:
            if line.startswith('!'):
                energy = float(line.split()[-2])    # in Ry
        magmom = np.nan    # not implemented yet...
    except:
        energy = np.nan    # error
        magmom = np.nan    # error
        print '    Structure ID {0}, could not obtain energy from {1}'.format(current_id, rin.outfile)

    #---------- collect the last structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(work_path+rin.outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(work_path+rin.infile)
        lines_atom = qe_structure.extract_atomic_positions(work_path+rin.outfile)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(work_path+rin.infile)
        opt_struc = qe_structure.from_lines(lines_cell, lines_atom)

        #----- opt_qe-structure
        with open('./data/opt_qe-structure', 'a') as fstruc:
            fstruc.write('# ID {0:d}\n'.format(current_id))
        qe_structure.write(opt_struc, './data/opt_qe-structure', mode='a')

    except:
        opt_struc = 'error'

    #---------- opt_CIFS
    if not opt_struc == 'error':
        cif = CifWriter(opt_struc, symprec=rin.symtoleR)
        cif.write_file(work_path+'tmp.cif')
        #----- correct title (need to delete '_chemical_formula_sum')
        with open(work_path+'tmp.cif', 'r') as fcif:
            ciflines = fcif.readlines()
        ciflines[1] = 'data_ID_{}\n'.format(current_id)
        if ciflines[11][:21] == '_chemical_formula_sum':
            ciflines.pop(11)
        else:
            raise ValueError('ciflines[11] is not _chemical_formula_sum, have to fix bag')
        #----- cif --> opt_cifs
        with open('./data/opt_CIFS.cif', 'a') as foptcif:
            for line in ciflines:
                foptcif.write(line)

    #---------- mv xxxxx fin_xxxxx
    qe_files = [rin.infile, 'tmp.cif', rin.outfile]
    for f in qe_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    #---------- clean stat file
    os.remove(work_path+'stat_job')

    #---------- return
    return opt_struc, energy, magmom, check_opt
