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
    return '1.4.0b10'


def set_logger(noprint=False, debug=False, logfile=None, errfile=None, debugfile=None):
    # ---------- level and formatter
    logger.setLevel(DEBUG)
    fmt = Formatter("[%(asctime)s][%(module)s][%(levelname)s] %(message)s")

    # ---------- stream handler for log
    if not noprint:
        shandler = StreamHandler(stream=sys.stdout)
        shandler.setFormatter(fmt)
        if debug:
            shandler.setLevel(DEBUG)
        else:
            shandler.setLevel(INFO)
        logger.addHandler(shandler)

    # ---------- file handler for log (level == INFO)
    if logfile is not None:
        fhandler = FileHandler(logfile)
        fhandler.setFormatter(fmt)
        fhandler.addFilter(lambda record: record.levelno == INFO)
        logger.addHandler(fhandler)

    # ---------- file handler for error (WARNING <= level)
    if errfile is not None:
        ehandler = FileHandler(errfile)
        ehandler.setFormatter(fmt)
        ehandler.setLevel(WARNING)
        logger.addHandler(ehandler)

    # ---------- file handler for debug (level == DEBUG)
    if debug:
        if debugfile is not None:
            dhandler = FileHandler(debugfile)
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


def backup_cryspy():

    # ---------- print
    logger.info('Backup data')

    # ---------- make directory
    dname = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = 'backup/' + dname + '/'
    os.makedirs(dst, exist_ok=True)

    # ---------- file/directory list
    flist = ['cryspy.in', 'cryspy.stat', 'debug_cryspy', 'err_cryspy', 'log_cryspy', 'cryspy_interactive.ipynb']
    dlist = ['calc_in', 'data']

    # ---------- backup
    for f in flist:
        if os.path.isfile(f):
            shutil.copy2(f, dst)
    for d in dlist:
        if os.path.isdir(d):
            shutil.copytree(d, dst + d)


def clean_cryspy(skip_yes=False):
    # ---------- yes/no
    if not skip_yes:
        while True:
            choice = input("Are you sure you want to clean the data? 'yes' or 'no' [y/n]: ").lower()
            if choice in ['y', 'yes']:
                break
            elif choice in ['n', 'no']:
                return

    # ---------- file/directory list
    fdlist = ['cryspy.stat', 'debug_cryspy', 'err_cryspy', 'log_cryspy', 'lock_cryspy',
              'data', 'work', 'tmp_gen_struc']

    # ---------- clean
    dname = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = 'trash/' + dname + '/'
    os.makedirs(dst, exist_ok=True)
    for fd in fdlist:
        if os.path.exists(fd):
            shutil.move(fd, dst)
