'''
Collect results in soiap
'''

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
                if ('QMD%frc converged.' in lines[i-2]
                        and 'QMD%strs converged.' in lines[i-1]):
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
        energy = energy/float(rin.natot)    # eV/cell --> eV/atom
    except:
        energy = np.nan    # error
        print('    Structure ID {0}, could not obtain energy from {1}'.format(
            current_id, rin.soiap_outfile))

    # ---------- collect the last structure
    try:
        opt_struc = soiap_structure.from_file(work_path+'log.struc')
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
