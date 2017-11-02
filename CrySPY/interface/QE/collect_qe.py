#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np
from pymatgen.core.units import Energy

from . import structure as qe_structure
from ...IO import read_input as rin


def collect_qe(current_id, work_path):
    # ---------- check optimization in previous stage
    try:
        with open(work_path+rin.qe_outfile, 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'End final coordinates' in line:
                check_opt = 'done'
    except:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    try:
        with open(work_path+rin.qe_outfile, 'r') as fpout:
            lines = fpout.readlines()
        energy = np.nan
        for line in reversed(lines):
            if line.startswith('!'):
                energy = float(line.split()[-2])    # in Ry
                energy = float(Energy(energy, 'Ry').to('eV'))    # Ry --> eV
                break
        magmom = np.nan    # not implemented yet...
    except:
        energy = np.nan    # error
        magmom = np.nan    # error
        print('    Structure ID {0}, could not obtain energy from {1}'.format(current_id, rin.qe_outfile))

    # ---------- collect the last structure
    try:
        lines_cell = qe_structure.extract_cell_parameters(work_path+rin.qe_outfile)
        if lines_cell is None:
            lines_cell = qe_structure.extract_cell_parameters(work_path+rin.qe_infile)
        lines_atom = qe_structure.extract_atomic_positions(work_path+rin.qe_outfile)
        if lines_atom is None:
            lines_atom = qe_structure.extract_atomic_positions(work_path+rin.qe_infile)
        opt_struc = qe_structure.from_lines(lines_cell, lines_atom)

        # ------ opt_qe-structure
        with open('./data/opt_qe-structure', 'a') as fstruc:
            fstruc.write('# ID {0:d}\n'.format(current_id))
        qe_structure.write(opt_struc, './data/opt_qe-structure', mode='a')

    except:
        opt_struc = None

    # ---------- mv xxxxx fin_xxxxx
    qe_files = [rin.qe_infile, rin.qe_outfile]
    for f in qe_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    # ---------- clean stat file
    os.remove(work_path+'stat_job')

    # ---------- return
    return opt_struc, energy, magmom, check_opt
