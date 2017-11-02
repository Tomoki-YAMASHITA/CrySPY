#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np
from pymatgen.core.units import Energy

from . import structure as soiap_structure
from ...IO import read_input as rin


def collect_soiap(current_id, work_path):
    # ---------- check optimization in current stage
    try:
        with open(work_path+rin.soiap_outfile, 'r') as fout:
            lines = fout.readlines()
        check_opt = 'not_yet'
        for i, line in enumerate(lines):
            if '*** QMD%loopc' in line:
                if 'QMD%frc converged.' in lines[i-2] and 'QMD%strs converged.' in lines[i-1]:
                    check_opt = 'done'
                break
    except:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    magmom = np.nan    # magnetic moment is not calculated
    try:
        with open(work_path+'log.tote') as f:
            lines = f.readlines()
        energy = float(lines[-1].split()[2])    # in Hartree
        energy = float(Energy(energy, 'Ha').to('eV'))    # Hartree --> eV
    except:
        energy = np.nan    # error
        print('    Structure ID {0}, could not obtain energy from {1}'.format(current_id, rin.soiap_outfile))

    # ---------- collect the last structure
    try:
        opt_struc = soiap_structure.from_file(work_path+'log.struc')
    except:
        opt_struc = None

    # ---------- mv xxxxx fin_xxxxx
    soiap_files = [rin.soiap_infile, rin.soiap_outfile, rin.soiap_cif,
                   'log.struc', 'log.tote', 'log.frc', 'log.strs']
    for f in soiap_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    # ---------- clean stat file
    os.remove(work_path+'stat_job')

    # ---------- return
    return opt_struc, energy, magmom, check_opt
