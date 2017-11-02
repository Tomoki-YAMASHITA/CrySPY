#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np
from pymatgen import Structure


def collect_vasp(current_id, work_path):
    # ---------- check optimization in previous stage
    try:
        with open(work_path+'prev_OUTCAR', 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'reached required accuracy' in line:
                check_opt = 'done'
    except:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom
    try:
        with open(work_path+'OSZICAR',  'r') as foszi:
            oszi = foszi.readlines()
        if 'F=' in oszi[-1]:
            energy = float(oszi[-1].split()[2])    # free energy (eV)
            if 'mag=' in oszi[-1]:
                magmom = float(oszi[-1].split()[-1])    # total magnetic moment
            else:
                magmom = np.nan
        else:
            energy = np.nan    # error
            magmom = np.nan    # error
            print('    Structure ID {0}, could not obtain energy from OSZICAR'.format(current_id))
    except:
        energy = np.nan    # error
        magmom = np.nan    # error
        print('    Structure ID {0}, could not obtain energy from OSZICAR'.format(current_id))

    # ---------- collect CONTCAR
    try:
        opt_struc = Structure.from_file(work_path+'CONTCAR')
    except:
        opt_struc = None

    # ---------- mv xxxxx fin_xxxxx
    vasp_files = ['POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR', 'WAVECAR', 'CHGCAR']
    for f in vasp_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    # ---------- remove STOPCAR
    if os.path.isfile(work_path+'STOPCAR'):
        os.remove(work_path+'STOPCAR')

    # ---------- clean stat file
    os.remove(work_path+'stat_job')

    # ---------- return
    return opt_struc, energy, magmom, check_opt
