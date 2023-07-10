'''
Calculate Fingerprint
'''

from logging import getLogger
import os
import subprocess

import numpy as np


logger = getLogger('cryspy')

class Calc_FP:
    '''
    calculate fingerprint using cal_fingerprint program

    # ---------- args
    struc_data (dict or list): structure data
        You may include None in struc_data
        if type of struc_data is list,
            struc_data is converted into dict type
            as {0: struc_data_0, 1: struc_data_1, ...}

    # ---------- instance methods
    self.calc(self)

    # ---------- comment
    descriptor data in self.descriptors
    '''

    def __init__(self, struc_data, fp_rmin=0.5, fp_rmax=5.0, fp_npoints=10,
                 fp_sigma=1.0, fppath='./cal_fingerprint'):
        # ---------- check args
        # ------ struc_data
        if type(struc_data) is dict:
            pass
        elif type(struc_data) is list:
            # -- convert to dict
            struc_data = {i: struc_data[i] for i in range(len(struc_data))}
        else:
            logger.error('Type of struc_data is wrong')
            raise SystemExit(1)
        self.struc_data = struc_data
        # ------ fp
        self.fp_rmin = fp_rmin
        self.fp_rmax = fp_rmax
        self.fp_sigma = fp_sigma
        self.fp_npoints = fp_npoints
        self.fppath = os.path.abspath(fppath)

    def calc(self):
        '''
        calculate fingerprint

        # ---------- return
        self.descriptors (dict): descriptor data
        '''
        # ---------- cd temporary directory
        if not os.path.isdir('tmp_calc_FP'):
            os.mkdir('tmp_calc_FP')
        os.chdir('tmp_calc_FP')
        # ---------- calc fingerprint
        self.descriptors = {}
        for cid, struc in self.struc_data.items():
            # ------ output POSCAR
            struc.to(fmt='poscar', filename='POSCAR')
            if not os.path.isfile('POSCAR'):
                logger.error('No POSCAR file')
                raise SystemExit(1)
            # ------ run cal_fingerprint
            with open('log_fingerprint', 'w') as logf:
                subprocess.call([self.fppath, 'POSCAR',
                                 '-rmin', '{}'.format(self.fp_rmin),
                                 '-rmax', '{}'.format(self.fp_rmax),
                                 '-npoints', '{}'.format(self.fp_npoints),
                                 '-sigma', '{}'.format(self.fp_sigma)],
                                stdout=logf, stderr=logf)
            fp = np.loadtxt('feature_ffpf.dat')
            # ------ mv xxx --> fin_xxx
            os.rename('POSCAR', 'fin_POSCAR')
            os.rename('feature_ffpf.dat', 'fin_feature_ffpf.dat')
            # ------ fp --> descriptors
            self.descriptors[cid] = fp
        # ---------- go back to ..
        os.chdir('../')
