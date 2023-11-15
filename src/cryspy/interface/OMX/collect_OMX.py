'''
Collect results in OpenMX
written  by  H. Sawahata 2020/03/09
info at hikaruri.jp
'''

from logging import getLogger

import numpy as np
import re
from pymatgen.core.units import Energy

from . import structure as OMX_structure
from ...IO import read_input as rin


logger = getLogger('cryspy')

def collect_OMX(current_id, work_path, nat):
    # ---------- check optimization in previous stage (done)
    # If *.out file exists, the calculation is done.
    try:
        with open(work_path+rin.OMX_outfile, 'r') as fpout:
            lines = fpout.readlines()
        check_opt = 'done'
    except Exception as e:
        logger.warning(str(e.args[0]))
        check_opt = 'no_file'

    # ---------- obtain energy and magmom (done)
    natot = sum(nat)    # do not use rin.natot here for EA-vc
    try:
        with open(work_path+rin.OMX_outfile, 'r') as fpout:
            lines = fpout.readlines()
        energy = np.nan
        for line in reversed(lines):
            if re.search('Utot\.', line):
                energy = float(re.search(r"-\d.+", line).group())
                energy = float(Energy(energy, 'Ha').to('eV'))
                energy = energy / float(natot)
                break
        magmom = np.nan # implemented (2020/10/03)
        for line in lines:
            if line.find("muB") >= 0:
                muB = line.split()
                magmom = float(muB[4])
                break
    except Exception as e:
        energy = np.nan    # error
        magmom = np.nan    # error
        logger.warning(str(e.args[0]) + f':    Structure ID {current_id},'
                       f' could not obtain energy from {rin.OMX_outfile}')

    # ---------- collect the last structure (yet)
    try:
        lines_cell = OMX_structure.extract_cell_parameters_from_outfile(
            work_path+rin.OMX_outfile)
        if lines_cell is None:
            lines_cell = OMX_structure.extract_cell_parameters_from_infile(
            work_path+rin.OMX_infile)
        lines_atom = OMX_structure.extract_atomic_positions_from_outfile(
            work_path+rin.OMX_outfile, nat)
        if lines_atom is None:
            lines_atom = OMX_structure.extract_atomic_positions_from_infile(
                work_path+rin.OMX_infile, nat)
        opt_struc = OMX_structure.from_lines(lines_cell, lines_atom)
        # ------ opt_OMX-structure
        with open('./data/opt_OMX-structure', 'a') as fstruc:
            fstruc.write('# ID {0:d}\n'.format(current_id))
        OMX_structure.write(opt_struc, './data/opt_OMX-structure', mode='a')
    except Exception as e:
        logger.warning(str(e.args[0]))
        opt_struc = None

    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan

    # ---------- return
    return opt_struc, energy, magmom, check_opt
