'''
Utility for CrySPY
'''

from datetime import datetime
import os
import shutil


# ---------- parameters
bohr2ang = 0.529177210903
hrt2ev = 27.211386245988
ry2ev = 13.605693122994
kbar2ev_ang3 = 0.0006241509073


# ---------- functions
def get_version():
    return 'CrySPY 0.11.0'


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


def backup_cryspy():
    # ---------- make directory
    dst = 'backup'
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.mkdir(dst)

    # ---------- file/directory list
    flist = ['cryspy.in', 'cryspy.out', 'cryspy.stat', 'err', 'log']
    dlist = ['calc_in', 'data', 'ext']

    # ---------- backup
    for f in flist:
        if os.path.isfile(f):
            shutil.copy2(f, dst)
    for d in dlist:
        if os.path.isdir(d):
            shutil.copytree(d, dst + '/' + d)
