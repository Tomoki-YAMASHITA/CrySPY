#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np

from .. import calc_dscrpt
from ..IO import read_input as rin


def calc_X(init_struc_data):
    if rin.dscrpt == 'FP':
        print('Calculate descriptors: FingerPrint')

        #---------- check cal_fingerprint executable file
        fppath = os.path.dirname(os.path.abspath(__file__)) + '/../f-fingerprint/cal_fingerprint'
        if not os.path.isfile(fppath):
            raise IOError('There is no cal_fingerprint program in CrySPY/f-fingerprint/cal_fingerprint')

        #---------- calc descriptors
        descriptors = calc_dscrpt.FP.calc_FP.calc_X(init_struc_data, fppath,
                                                    rin.fp_rmin, rin.fp_rmax,
                                                    rin.fp_npoints, rin.fp_sigma)
    else:
        raise ValueError('Now FP only')

    return descriptors


def append_X(init_struc_data, next_BO_id, non_error_id, descriptors):
    if rin.dscrpt == 'FP':
        print('Append descriptors: FingerPrint')

        #---------- check cal_fingerprint executable file
        fppath = os.path.dirname(os.path.abspath(__file__)) + '/../f-fingerprint/cal_fingerprint'
        if not os.path.isfile(fppath):
            raise IOError('There is no cal_fingerprint program in CrySY/f-fingerprint/cal_fingerprint')

        #---------- append descriptor
        tmp_dscrpt = calc_dscrpt.FP.calc_FP.calc_X(init_struc_data[next_BO_id:], fppath,
                                                   rin.fp_rmin, rin.fp_rmax,
                                                   rin.fp_npoints, rin.fp_sigma)
        descriptors = np.vstack((descriptors, tmp_dscrpt))
        non_error_id = np.r_[non_error_id, np.arange(next_BO_id, len(init_struc_data))]

    else:
        raise ValueError('Now FP only')

    return non_error_id, descriptors
