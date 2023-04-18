'''
Calculation files in ASE
'''

import os

from ...IO import read_input as rin


def check_input_ase():
    # ---------- prepare rin.jobfile, POTCAR, INCAR
    calc_inputs = [rin.jobfile, rin.ase_python]

    # ----- check required files
    for f in calc_inputs:
        if f == rin.ase_python:
            finfiles = [rin.ase_python + '_{}'.format(i) for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    raise IOError('Could not find ./calc_in/'+ff)
        else:
            if not os.path.isfile('./calc_in/'+f):
                raise IOError('Could not find ./calc_in/'+f)