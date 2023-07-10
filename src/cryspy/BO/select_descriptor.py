'''
Select descriptors
    - FP: FingerPrint
'''

from logging import getLogger

from ..util.utility import check_fppath
from ..calc_dscrpt.FP.calc_FP import Calc_FP
from ..IO import read_input as rin


logger = getLogger('cryspy')

def select_descriptor(struc_data):
    # ---------- fingerprint
    if rin.dscrpt == 'FP':
        logger.info('Calculate descriptors: FingerPrint')
        # ------ check cal_fingerprint executable file
        fppath = check_fppath(rin.fppath)
        # ------ calc fingerprint
        cfp = Calc_FP(struc_data, rin.fp_rmin, rin.fp_rmax,
                      rin.fp_npoints, rin.fp_sigma, fppath)
        cfp.calc()
        return cfp.descriptors
    else:
        logger.error('Now FP only')
        raise SystemExit(1)
