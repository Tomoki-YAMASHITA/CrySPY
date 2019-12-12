#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from .. import utility
from ..calc_dscrpt.FP.calc_FP import Calc_FP
from ..IO import read_input as rin


def select_descriptor(struc_data):
    # ---------- fingerprint
    if rin.dscrpt == 'FP':
        print('Calculate descriptors: FingerPrint')
        # ------ check cal_fingerprint executable file
        fppath = utility.check_fppath()
        # ------ calc fingerprint
        cfp = Calc_FP(struc_data, rin.fp_rmin, rin.fp_rmax,
                      rin.fp_npoints, rin.fp_sigma, fppath)
        cfp.calc()
        return cfp.descriptors
    else:
        raise NotImplementedError('Now FP only')
