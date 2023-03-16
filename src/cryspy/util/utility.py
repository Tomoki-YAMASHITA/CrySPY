'''
Utility for CrySPY
'''

from datetime import datetime
import os
import shutil
import subprocess


# ---------- functions
def get_version():
    return '1.0.0'


def get_date():
    return datetime.now().strftime("%Y/%m/%d %H:%M:%S")


def check_fwpath(fwpath):
    if fwpath is None:
        # ---------- check if find_wy is in your path
        sr = subprocess.run(['which', 'find_wy'], capture_output=True, text=True)
        fwpath = sr.stdout.strip()    # to delete \n
        if fwpath == '':
            raise IOError('There is no find_wy program in your path')
    else:
        # ---------- check fwpath written in cryspy.in
        if not os.path.isfile(fwpath):
            raise IOError('There is no find_wy program in {}'.format(fwpath))
    return fwpath


def check_fppath(fppath):
    if fppath is None:
        # ---------- check if find_wy is in your path
        sr = subprocess.run(['which', 'cal_fingerprint'], capture_output=True, text=True)
        fppath = sr.stdout.strip()    # to delete \n
        if fppath == '':
            raise IOError('There is no cal_fingerprint program in your path')
    else:
        # ---------- check fppath written in cryspy.in
        if not os.path.isfile(fppath):
            raise IOError('There is no cal_fingerprint program in {}'.format(fppath))
    return fppath


def backup_cryspy():
    # ---------- make directory
    dst = 'backup'
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.mkdir(dst)

    # ---------- file/directory list
    flist = ['cryspy.in', 'cryspy.stat', 'err_cryspy', 'log_cryspy']
    dlist = ['calc_in', 'data', 'ext']

    # ---------- backup
    for f in flist:
        if os.path.isfile(f):
            shutil.copy2(f, dst)
    for d in dlist:
        if os.path.isdir(d):
            shutil.copytree(d, dst + '/' + d)

    # ---------- print
    print('\nBackup data')


def clean_cryspy():
    # ---------- yes/no
    while True:
        choice = input("Are you sure you want to clean the data? 'yes' or 'no' [y/n]: ").lower()
        if choice in ['y', 'yes']:
            break
        elif choice in ['n', 'no']:
            return

    # ---------- file/directory list
    fdlist = ['cryspy.stat', 'err_cryspy', 'log_cryspy', 'lock_cryspy',
              'data', 'work', 'ext', 'tmp_calc_FP', 'tmp_gen_struc']

    # ---------- clean
    dname = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = 'trash/' + dname + '/'
    os.makedirs(dst, exist_ok=True)
    for fd in fdlist:
        if os.path.exists(fd):
            shutil.move(fd, dst)
