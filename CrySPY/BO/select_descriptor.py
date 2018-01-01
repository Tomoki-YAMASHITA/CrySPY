#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os

import numpy as np

from ..calc_dscrpt.FP import calc_FP
from ..IO import read_input as rin


def calc_X(init_struc_data):
    if rin.dscrpt == 'FP':
        print('Calculate descriptors: FingerPrint')

        # ---------- check cal_fingerprint executable file
        fppath = os.path.dirname(os.path.abspath(__file__)) + '/../f-fingerprint/cal_fingerprint'
        if not os.path.isfile(fppath):
            raise IOError('There is no cal_fingerprint program in CrySPY/f-fingerprint/cal_fingerprint')

        # ---------- calc descriptors
        strucs_list = [init_struc_data[i] for i in range(len(init_struc_data))]    # dict --> list
        descriptors = calc_FP.calc_X(strucs_list, fppath,
                                     rin.fp_rmin, rin.fp_rmax,
                                     rin.fp_npoints, rin.fp_sigma)
    else:
        raise ValueError('Now FP only')
    return descriptors


def append_X(init_struc_data, next_BO_id, non_error_id, descriptors):
    if rin.dscrpt == 'FP':
        print('Append descriptors: FingerPrint')

        # ---------- check cal_fingerprint executable file
        fppath = os.path.dirname(os.path.abspath(__file__)) + '/../f-fingerprint/cal_fingerprint'
        if not os.path.isfile(fppath):
            raise IOError('There is no cal_fingerprint program in CrySY/f-fingerprint/cal_fingerprint')

        # ---------- append descriptor
        strucs_list = [init_struc_data[i] for i in range(next_BO_id, len(init_struc_data))]    # dict --> list
        tmp_dscrpt = calc_FP.calc_X(strucs_list, fppath,
                                    rin.fp_rmin, rin.fp_rmax,
                                    rin.fp_npoints, rin.fp_sigma)
        descriptors = np.vstack((descriptors, tmp_dscrpt))
        non_error_id = np.r_[non_error_id, np.arange(next_BO_id, len(init_struc_data))]
    else:
        raise ValueError('Now FP only')
    return non_error_id, descriptors
