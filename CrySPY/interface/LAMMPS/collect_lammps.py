#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np

from . import structure as lammps_structure
from ...IO import read_input as rin


def collect_lammps(current_id, work_path):
    # ---------- check optimization in current stage & obtain energy
    energy = np.nan
    check_opt = 'not_yet'
    try:
        with open(work_path+rin.lammps_outfile, 'r') as fout:
            lines = fout.readlines()
        for i, line in enumerate(lines):
            if 'ERROR:' in line:
                print('    ' + line.rstrip())
                energy = np.nan
                check_opt = 'not_yet'
                break
            elif 'Minimization stats:' in line:    # Provisional judgment condition
                energy = float(lines[i+3].split()[2])    # in eV (if units is metal)
                energy = energy/float(rin.natot)    # eV/cell --> eV/atom
                check_opt = 'done'
    except:
        print('    Structure ID {0}, could not obtain energy from {1}'.format(current_id, rin.lammps_outfile))
        energy = np.nan    # error
        check_opt = 'no_file'

    # ---------- obtain magmom
    magmom = np.nan    # magnetic moment is not calculated

    # ---------- collect the last structure
    try:
        opt_struc = lammps_structure.from_file(work_path+'log.struc')
    except:
        opt_struc = None
    # ---------- mv xxxxx fin_xxxxx
    if rin.lammps_potential is None:
        lammps_files = [rin.lammps_infile, rin.lammps_outfile,
                        rin.lammps_data, 'log.struc']
    else:
        lammps_files = [rin.lammps_infile, rin.lammps_outfile,
                        rin.lammps_potential, rin.lammps_data, 'log.struc']
    for f in lammps_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    # ---------- clean stat file
    os.remove(work_path+'stat_job')

    # ---------- return
    return opt_struc, energy, magmom, check_opt
