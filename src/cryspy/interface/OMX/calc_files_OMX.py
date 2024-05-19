'''
Calculation files in OpenMX
written by H. Sawahata 2020/03/09
info at hikaruri.jp
'''

from logging import getLogger
import os


logger = getLogger('cryspy')


def check_input_OMX(rin):
    # ---------- prepare rin.jobfile, rin.OMX_infile
    calc_inputs = [rin.jobfile, rin.OMX_infile]

    # ------ check required files
    for f in calc_inputs:
        if f == rin.OMX_infile:
            finfiles = [rin.OMX_infile + f'_{i}' for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/' + ff):
                    logger.error('Could not find ./calc_in/' + ff)
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                raise SystemExit(1)
