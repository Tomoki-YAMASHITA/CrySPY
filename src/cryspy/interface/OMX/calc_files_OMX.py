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
            for i in range(1, rin.nstage+1):
                fname_candidates = [
                    f'{i}_{rin.OMX_infile}',
                    f'{rin.OMX_infile}_{i}',
                    f'{rin.OMX_infile}'
                ]
                found = False
                for fname in fname_candidates:
                    if os.path.isfile('./calc_in/' + fname):
                        found = True
                        break
                if not found:
                    logger.error('Could not find in ./calc_in/: ' + fname_candidates[0] + ' or ' + fname_candidates[-1])
                    os.remove('lock_cryspy')
                    raise SystemExit(1)
        else:
            if not os.path.isfile('./calc_in/' + f):
                logger.error('Could not find ./calc_in/' + f)
                os.remove('lock_cryspy')
                raise SystemExit(1)
