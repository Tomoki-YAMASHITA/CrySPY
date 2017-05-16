#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np
from pymatgen import Structure
from pymatgen.io.cif import CifWriter

from ...IO import read_input as rin


def collect_vasp(current_id, work_path):
    #---------- check optimization in previous stage
    try:
        with open(work_path+'prev_OUTCAR', 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'not_yet'
        for line in lines:
            if 'reached required accuracy' in line:
                check_opt = 'done'
    except:
        check_opt = 'no_file'

    #---------- obtain energy and magmom
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

    #---------- collect CONTCAR
    try:
        opt_struc = Structure.from_file(work_path+'CONTCAR')

#        #----- opt_POSCARS
#        with open(work_path+'CONTCAR', 'r') as ft:
#            lines = ft.readlines()
#        with open('./data/opt_POSCARS', 'a') as fopt:
#            for line in lines:
#                if len(line) == 2:    # only '\n'
#                    break
#                else:
#                    fopt.write(line)
    except:
        opt_struc = None

#    #---------- opt_CIFS
#    if is not None:
#        cif = CifWriter(opt_struc, symprec=rin.symtoleR)
#        cif.write_file(work_path+'tmp.cif')
#        #----- correct title (need to delete '_chemical_formula_sum')
#        with open(work_path+'tmp.cif', 'r') as fcif:
#            ciflines = fcif.readlines()
#        ciflines[1] = 'data_ID_{}\n'.format(current_id)
#        if ciflines[11][:21] == '_chemical_formula_sum':
#            ciflines.pop(11)
#        else:
#            raise ValueError('ciflines[11] is not _chemical_formula_sum, have to fix bag')
#        #----- cif --> opt_cifs
#        with open('./data/opt_CIFS.cif', 'a') as foptcif:
#            for line in ciflines:
#                foptcif.write(line)

    #---------- mv xxxxx fin_xxxxx
    vasp_files = ['POSCAR', 'CONTCAR', 'OUTCAR', 'OSZICAR', 'WAVECAR', 'CHGCAR']
    for f in vasp_files:
        if os.path.isfile(work_path+f):
            os.rename(work_path+f, work_path+'fin_'+f)

    #---------- clean stat file
    os.remove(work_path+'stat_job')

    #---------- return
    return opt_struc, energy, magmom, check_opt

