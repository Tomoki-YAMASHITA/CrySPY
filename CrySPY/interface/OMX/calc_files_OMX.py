'''
Calculation files in OpenMX
written by H. Sawahata 2020/03/09
info at hikaruri.jp
'''

import os

from ...IO import read_input as rin


def check_input_OMX():
    # ---------- prepare rin.jobfile, rin.OMX_infile
    calc_inputs = [rin.jobfile, rin.OMX_infile]

    # ------ check required files
    for f in calc_inputs:
        if f == rin.OMX_infile:
            finfiles = [rin.OMX_infile + '_{}'.format(i) for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/' + ff):
                    raise IOError('Could not find ./calc_in/' + ff)
        else:
            if not os.path.isfile('./calc_in/' + f):
                raise IOError('Could not find ./calc_in/' + f)
