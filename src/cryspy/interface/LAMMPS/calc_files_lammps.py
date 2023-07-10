'''
Calculation files in LAMMPS
'''

from logging import getLogger
import os

from ...IO import read_input as rin


logger = getLogger('cryspy')

def check_input_lammps():
    # ---------- prepare rin.jobfile, rin.lammps_potential, rin.rammps_infile
    if rin.lammps_potential is None:
        calc_inputs = [rin.jobfile, rin.lammps_infile]
    else:
        calc_inputs = [rin.jobfile, rin.lammps_infile] + rin.lammps_potential

    # ----- check required files
    for f in calc_inputs:
        if f == rin.lammps_infile:
            finfiles = [rin.lammps_infile + '_{}'.format(i) for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    logger.error('Could not find ./calc_in/'+ff)
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/'+f):
                logger.error('Could not find ./calc_in/'+f)
                raise SystemExit(1)