from logging import getLogger
import numpy as np

from . import structure as lammps_structure


logger = getLogger('cryspy')


def collect_lammps(rin, cid, work_path, nat):
    # ---------- natot
    natot = sum(nat)
    # ---------- check optimization in current stage & obtain energy
    energy = np.nan
    check_opt = 'not_yet'
    try:
        with open(work_path+rin.lammps_outfile, 'r') as fout:
            lines = fout.readlines()
        for i, line in enumerate(lines):
            if 'ERROR:' in line:
                logger.warning('    ' + line.rstrip())
                energy = np.nan
                check_opt = 'not_yet'
                break
            elif 'Minimization stats:' in line:
                energy = float(lines[i+3].split()[2])  # in eV (units is metal)
                energy = energy/float(natot)    # eV/cell --> eV/atom
                check_opt = 'done'
    except Exception as e:
        logger.warning(f'{e}:    Structure ID {cid},'
                       f'could not obtain energy from {rin.lammps_outfile}')
        energy = np.nan    # error
        check_opt = 'no_file'

    # ---------- obtain magmom
    magmom = np.nan    # magnetic moment is not calculated

    # ---------- collect the last structure
    try:
        opt_struc = lammps_structure.from_file(rin, work_path+'log.struc', nat)
    except Exception as e:
        logger.warning(e)
        opt_struc = None

    # ---------- return
    return opt_struc, energy, magmom, check_opt
