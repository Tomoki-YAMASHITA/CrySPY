'''
Utility for CrySPY
'''

from datetime import datetime
from logging import getLogger, StreamHandler, FileHandler, Formatter, DEBUG, INFO, WARNING
import os
import shutil
import subprocess
import sys


logger = getLogger('cryspy')

def get_version():
    return '1.2.3'


def set_logger(noprint=False, debug=False):
    # ---------- level and formatter
    logger.setLevel(DEBUG)
    fmt = Formatter("[%(asctime)s][%(module)s][%(levelname)s] %(message)s")

    # ---------- stream handler for log
    if not noprint:
        shandler = StreamHandler(stream=sys.stdout)
        shandler.setFormatter(fmt)
        shandler.setLevel(INFO)
        logger.addHandler(shandler)

    # ---------- file handler for log (level == INFO)
    fhandler = FileHandler('./log_cryspy')
    fhandler.setFormatter(fmt)
    fhandler.addFilter(lambda record: record.levelno == INFO)
    logger.addHandler(fhandler)

    # ---------- file handler for error (WARNING <= level)
    ehandler = FileHandler('./err_cryspy')
    ehandler.setFormatter(fmt)
    ehandler.setLevel(WARNING)
    logger.addHandler(ehandler)

    # ---------- file handler for debug (level == DEBUG)
    if debug:
        dhandler = FileHandler('./debug_cryspy')
        dhandler.setFormatter(fmt)
        dhandler.addFilter(lambda record: record.levelno == DEBUG)
        logger.addHandler(dhandler)


def check_fwpath(fwpath):
    if fwpath is None:
        # ---------- check if find_wy is in your path
        sr = subprocess.run(['which', 'find_wy'], capture_output=True, text=True)
        fwpath = sr.stdout.strip()    # to delete \n
        if fwpath == '':
            logger.error('There is no find_wy program in your path')
            raise SystemExit(1)
    else:
        # ---------- check fwpath written in cryspy.in
        if not os.path.isfile(fwpath):
            logger.error(f'There is no find_wy program in {fwpath}')
            raise SystemExit(1)
    return fwpath


def check_fppath(fppath):
    if fppath is None:
        # ---------- check if find_wy is in your path
        sr = subprocess.run(['which', 'cal_fingerprint'], capture_output=True, text=True)
        fppath = sr.stdout.strip()    # to delete \n
        if fppath == '':
            logger.error('There is no cal_fingerprint program in your path')
            raise SystemExit(1)
    else:
        # ---------- check fppath written in cryspy.in
        if not os.path.isfile(fppath):
            logger.error(f'There is no cal_fingerprint program in {fppath}')
            raise SystemExit(1)
    return fppath


def backup_cryspy():
    # ---------- make directory
    dname = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = 'backup/' + dname + '/'
    os.makedirs(dst, exist_ok=True)

    # ---------- file/directory list
    flist = ['cryspy.in', 'cryspy.stat', 'err_cryspy', 'log_cryspy']
    dlist = ['calc_in', 'data', 'ext']

    # ---------- backup
    for f in flist:
        if os.path.isfile(f):
            shutil.copy2(f, dst)
    for d in dlist:
        if os.path.isdir(d):
            shutil.copytree(d, dst + d)

    # ---------- print
    logger.info('Backup data')


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
