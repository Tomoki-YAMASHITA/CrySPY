'''
Collect results in OpenMX
written  by  H. Sawahata 2020/03/09
info at hikaruri.jp
'''

import numpy as np
import re
from pymatgen.core.units import Energy

from . import structure as OMX_structure
from ...IO import read_input as rin


def collect_OMX(current_id, work_path):
    # ---------- check optimization in previous stage (done)
    # If *.out file exists, the calculation is done.
    try:
        with open(work_path+rin.OMX_outfile, 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'done'
    except:
        check_opt = 'no_file'

    # ---------- obtain energy and magmom (done)
    try:
        with open(work_path+rin.OMX_outfile, 'r') as fpout:
            lines = fpout.readlines()
        energy = np.nan
        for line in reversed(lines):
            if re.search('Utot\.', line):
                energy = float(re.search(r"-\d.+", line).group())
                energy = float(Energy(energy, 'Ha').to('eV'))
                energy = energy / float(rin.natot)
                break	
        magmom = np.nan # not implement (2020/03/05)
    except:
        energy = np.nan    # error
        magmom = np.nan    # error
        print(' Structure ID {0}, could not obtain energy from {1}'.format(
            current_id, rin.OMX_outfile))

    # ---------- collect the last structure (yet)
    try:
        lines_cell = OMX_structure.extract_cell_parameters_from_outfile(
            work_path+rin.OMX_outfile)
        if lines_cell is None:
            lines_cell = OMX_structure.extract_cell_parameters_from_infile(
            work_path+rin.OMX_infile)
        lines_atom = OMX_structure.extract_atomic_positions_from_outfile(
            work_path+rin.OMX_outfile)
        if lines_atom is None:
            lines_atom = OMX_structure.extract_atomic_positions_from_infile(
                work_path+rin.OMX_infile)
        opt_struc = OMX_structure.from_lines(lines_cell, lines_atom)
        # ------ opt_OMX-structure
        with open('./data/opt_OMX-structure', 'a') as fstruc:
            fstruc.write('# ID {0:d}\n'.format(current_id))
        OMX_structure.write(opt_struc, './data/opt_OMX-structure', mode='a')
    except:
        opt_struc = None

    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan

    # ---------- return
    return opt_struc, energy, magmom, check_opt
