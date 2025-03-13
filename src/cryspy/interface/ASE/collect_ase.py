from logging import getLogger
import warnings

import numpy as np
from pymatgen.core import Structure


logger = getLogger('cryspy')


def collect_ase(cid, work_path, nat):
    # ---------- natot
    natot = sum(nat)
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
        logger.warning(f'{e}:    Structure ID {cid}, could not obtain energy from log.tote')
    # ---------- collect CONTCAR
    try:
        with warnings.catch_warnings():
            # to avoid BadPoscarWarning message from Pymatgen
            warnings.simplefilter('ignore')
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