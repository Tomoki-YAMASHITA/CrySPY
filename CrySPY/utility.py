#!/usr/bin/env python
# -*- coding: utf-8 -*-


from datetime import datetime
import os


def get_version():
    return 'CrySPY 0.6.1'


def get_date():
    return datetime.now().strftime("%Y/%m/%d %H:%M:%S")


def get_init_pos_path():
    return '../data/init_POSCARS'


def check_fwpath():
    # ---------- check find_wy executable file
    fwpath = os.path.dirname(os.path.abspath(__file__)) + '/find_wy/find_wy'
    if not os.path.isfile(fwpath):
        raise IOError('There is no find_wy program in {}'.format(fwpath))
    return fwpath
