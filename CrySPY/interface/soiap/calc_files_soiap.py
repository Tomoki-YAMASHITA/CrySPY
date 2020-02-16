'''
Calculation files in soiap
'''

import os

from ...IO import read_input as rin


def check_input_soiap():
    # ---------- prepare rin.jobfile, rin.soiap_infile
    calc_inputs = [rin.jobfile, rin.soiap_infile]

    # ----- check required files
    for f in calc_inputs:
        if f == rin.soiap_infile:
            finfiles = [rin.soiap_infile + '_{}'.format(i)
                        for i in range(1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    raise IOError('Could not find ./calc_in/'+ff)
        else:
            if not os.path.isfile('./calc_in/'+f):
                raise IOError('Could not find ./calc_in/'+f)
