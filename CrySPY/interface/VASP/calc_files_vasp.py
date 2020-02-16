'''
Calculation files in VASP
'''

import os

from ...IO import read_input as rin


def check_input_vasp():
    # ---------- prepare rin.jobfile, POTCAR, INCAR
    calc_inputs = [rin.jobfile, 'POTCAR', 'INCAR']

    # ------ check required files
    for f in calc_inputs:
        if f == 'INCAR':
            fincars = ['INCAR_{}'.format(i) for i in range(1, rin.nstage+1)]
            for ff in fincars:
                if not os.path.isfile('./calc_in/' + ff):
                    raise IOError('Could not find ./calc_in/' + ff)
        else:
            if not os.path.isfile('./calc_in/' + f):
                raise IOError('Could not find ./calc_in/' + f)
