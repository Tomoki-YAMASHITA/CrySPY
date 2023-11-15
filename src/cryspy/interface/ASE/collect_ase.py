'''
Collect results in ASE
'''

from logging import getLogger

import numpy as np
from pymatgen.core import Structure

from ...IO import read_input as rin


logger = getLogger('cryspy')

def collect_ase(current_id, work_path, nat):
    # ---------- natot
    natot = sum(nat)    # do not use rin.natot here for EA-vc
    # ---------- etc
    magmom = np.nan
    check_opt = 'no_file'    # always no_file in ASE for now
    # ---------- collect energy
    try:
        with open(work_path+'log.tote') as f:
            lines = f.readlines()
        energy = float(lines[-1].split()[0])    # in eV/cell
        energy = energy/float(natot)    # eV/cell --> eV/atom
    except Exception as e:
        energy = np.nan    # error
        logger.warning(str(e.args[0]) + f':    Structure ID {current_id}, could not obtain energy from log.tote')
    # ---------- collect CONTCAR
    try:
        opt_struc = Structure.from_file(work_path+'CONTCAR')
    except Exception:
        opt_struc = None
    # ---------- check
    if np.isnan(energy):
        opt_struc = None
    if opt_struc is None:
        energy = np.nan
        magmom = np.nan
    # ---------- return
    return opt_struc, energy, magmom, check_opt