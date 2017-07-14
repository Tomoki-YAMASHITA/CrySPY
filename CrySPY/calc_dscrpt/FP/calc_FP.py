#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess

import numpy as np


def calc_X(strucs_list, fppath='./cal_fingerprint', fp_rmin=0.5, fp_rmax=5.0, fp_npoints=50, fp_sigma=0.2):
    #---------- cd calc_descriptor
    if not os.path.isdir('calc_FP'):
        os.mkdir('calc_FP')
    os.chdir('calc_FP')

    #---------- calc fingerprint
    descriptors = None
    for struc in strucs_list:
        #------ output POSCAR
        struc.to(fmt='poscar', filename='POSCAR')
        fp = calc_fingerprint(fppath, fp_rmin, fp_rmax, fp_npoints, fp_sigma)
        if descriptors is None:
            descriptors = fp
        else:
            descriptors = np.vstack((descriptors, fp))

    #---------- go back to ..
    os.chdir('../')

    return descriptors


def calc_fingerprint(fppath='./cal_fingerprint', fp_rmin=0.5, fp_rmax=5.0, fp_npoints=50, fp_sigma=0.2):
    #---------- calc fingerprint
    if not os.path.isfile('POSCAR'):
        raise IOError('No POSCAR file')
    with open('sublog', 'w') as logf:
        subprocess.call([fppath, 'POSCAR', '-rmin', '{}'.format(fp_rmin), '-rmax', '{}'.format(fp_rmax),
                         '-npoints', '{}'.format(fp_npoints), '-sigma', '{}'.format(fp_sigma)],
                        stdout=logf, stderr=logf)
    fp = np.loadtxt('feature_ffpf.dat')

    #---------- mv xxx --> fin_xxx
    os.rename('POSCAR', 'fin_POSCAR')
    os.rename('feature_ffpf.dat', 'fin_feature_ffpf.dat')

    return fp
