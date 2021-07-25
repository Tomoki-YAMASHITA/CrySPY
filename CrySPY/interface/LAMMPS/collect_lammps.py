'''
Collect results in LAMMPS
'''

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
            elif 'Minimization stats:' in line:
                energy = float(lines[i+3].split()[2])  # in eV (units is metal)
                energy = energy/float(rin.natot)    # eV/cell --> eV/atom
                check_opt = 'done'
    except Exception as e:
        print(e)
        print('    Structure ID {0}, could not obtain energy from {1}'.format(
            current_id, rin.lammps_outfile))
        energy = np.nan    # error
        check_opt = 'no_file'

    # ---------- obtain magmom
    magmom = np.nan    # magnetic moment is not calculated

    # ---------- collect the last structure
    try:
        opt_struc = lammps_structure.from_file(work_path+'log.struc')
    except Exception as e:
        print(e)
        opt_struc = None

    # ---------- return
    return opt_struc, energy, magmom, check_opt
