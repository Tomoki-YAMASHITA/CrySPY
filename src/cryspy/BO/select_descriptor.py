from logging import getLogger

from ..calc_dscrpt.FP.calc_FP import calc_fp


logger = getLogger('cryspy')


def select_descriptor(rin, struc_data):
    # ---------- fingerprint
    if rin.dscrpt == 'FP':
        logger.info('Calculate descriptors: FingerPrint by Valle and Oganov')
        # ------ calc fingerprint
        descriptors = calc_fp(struc_data, rin.atype, rin.fp_rmax,
                              rin.fp_npoints, rin.fp_sigma)
        return descriptors
    else:
        logger.error('Now FP only')
        raise SystemExit(1)
