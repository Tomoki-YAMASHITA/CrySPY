'''
Utility for CrySPY
'''

from datetime import datetime
import os


# ---------- parameters
bohr2ang = 0.529177210903
hrt2ev = 27.211386245988
kbar2ev_ang3 = 0.0006241509073

# ---------- functions
def get_version():
    return 'CrySPY 0.9.0'


def get_date():
    return datetime.now().strftime("%Y/%m/%d %H:%M:%S")


def check_fwpath():
    # ---------- check find_wy executable file
    fwpath = os.path.dirname(os.path.abspath(__file__)) + '/find_wy/find_wy'
    if not os.path.isfile(fwpath):
        raise IOError('There is no find_wy program in {}'.format(fwpath))
    return fwpath


def check_fppath():
    # ---------- check cal_fingerprint executable file
    fppath = os.path.dirname(
        os.path.abspath(__file__)) + '/f-fingerprint/cal_fingerprint'
    if not os.path.isfile(fppath):
        raise IOError('There is no cal_fingerprint program in {}'.format(
            fppath))
    return fppath
