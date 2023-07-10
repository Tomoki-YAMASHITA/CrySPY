'''
Calculation files in Quantum ESPRESSO
'''

from logging import getLogger
import os

from ...IO import read_input as rin


logger = getLogger('cryspy')

def check_input_qe():
    # ---------- prepare rin.jobfile, rin.qe_infile
    calc_inputs = [rin.jobfile, rin.qe_infile]

    # ------ check required files
    for f in calc_inputs:
        if f == rin.qe_infile:
            finfiles = [rin.qe_infile + '_{}'.format(i) for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/' + ff):
                    logger.error('Could not find ./calc_in/' + ff)
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                raise SystemExit(1)