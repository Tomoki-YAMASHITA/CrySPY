#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import subprocess

from gen_cell import rndgen_lattice, calc_latvec, calc_cos
from with_spg.fw import fw_input, gen_wypos
from wo_spg.gen_coordinate import rndgen_coord
from ..struc_util import out_poscar


def rndgen_wo_spg(nstruc, natot, atype, nat, id_offset=0, minlen=4, maxlen=10, dangle=20, mindist=1.5,
                  maxcnt=200, symprec=0.001, init_pos_path=None):
    '''
    Randomly generate structures without space group information
    '''

    # ---------- initialize
    init_struc_data = {}
    spgnum = 0

    # ---------- cd gen_struc
    if not os.path.isdir('gen_struc'):
        os.mkdir('gen_struc')
    os.chdir('gen_struc')

    # ---------- generate structures
    while len(init_struc_data) < nstruc:
        spg_in, a, b, c, alpha, beta, gamma = rndgen_lattice(spgnum, minlen, maxlen, dangle)
        va, vb, vc = calc_latvec(a, b, c, alpha, beta, gamma)
        tmp_struc = rndgen_coord(natot, atype, nat, va, vb, vc, mindist, maxcnt)
        if tmp_struc is not None:    # success of generation
            # ------ check actual space group using pymatgen
            try:
                spg_sym, spg_num = tmp_struc.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # ------ register the structure in pymatgen format
            cID = len(init_struc_data) + id_offset
            init_struc_data[cID] = tmp_struc
            print('Structure ID {0:>8} was generated. Space group: {1:>3} --> {2:>3} {3}'.format(
                   cID, spg_in, spg_num, spg_sym))
            # ------ save poscar
            if init_pos_path is not None:
                out_poscar(tmp_struc, cID, init_pos_path)
    # ---------- go back to ..
    os.chdir('../')

    return init_struc_data


def rndgen_spg(nstruc, natot, atype, nat, spgnum='all', id_offset=0,
               minlen=4, maxlen=10, dangle=20, mindist=1.5,
               maxcnt=200, symprec=0.001,
               init_pos_path=None, fwpath='./find_wy'):
    '''
    Randomly generate structures with space group information
    '''

    # ---------- initialize
    init_struc_data = {}

    # ---------- cd gen_struc
    if not os.path.isdir('gen_struc'):
        os.mkdir('gen_struc')
    os.chdir('gen_struc')

    # ---------- generate structures
    while len(init_struc_data) < nstruc:
        spg_in, a, b, c, alpha, beta, gamma = rndgen_lattice(spgnum, minlen, maxlen, dangle)
        cosa, cosb, cosg = calc_cos(alpha, beta, gamma)
        fw_input(atype, nat, spg_in, a, b, c, cosa, cosb, cosg)

        # ------ loop for same fw_input
        cnt = 0
        while cnt <= maxcnt:
            # -- run find_wy
            with open('sublog', 'w') as f:
                subprocess.call([fwpath, 'input'], stdout=f, stderr=f)

            # -- generate a structure using POS_WY_SKEL_ALL.json
            if not os.path.isfile('POS_WY_SKEL_ALL.json'):
                wyflag = False
                break
            wyflag, tmp_struc = gen_wypos(atype, mindist, maxcnt)
            if wyflag is False:    # Failure
                os.remove('POS_WY_SKEL_ALL.json')
                cnt += 1
                continue
            else:    # Success
                rm_files()    # clean
                break         # break fw_input loop

        if wyflag is False:    # maximum trial or no POS_WY_SKEL_ALL.json file
            rm_files()    # clean
            continue      # to new fw_input

        # ------ check actual space group using pymatgen
        try:
            spg_sym, spg_num = tmp_struc.get_space_group_info(symprec=symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None

        # ------ register the structure in pymatgen format
        cID = len(init_struc_data) + id_offset
        init_struc_data[cID] = tmp_struc
        print('Structure ID {0:>8} was generated. Space group: {1:>3} --> {2:>3} {3}'.format(
               cID, spg_in, spg_num, spg_sym))

        # ------ save poscar
        if init_pos_path is not None:
            out_poscar(tmp_struc, cID, init_pos_path)

        # ------ clean
        rm_files()

    # ---------- go back to ..
    os.chdir('../')

    return init_struc_data


def rm_files(files=['input', 'POS_WY_SKEL_ALL.json']):
    for rfile in files:
        if os.path.isfile(rfile):
            os.remove(rfile)
