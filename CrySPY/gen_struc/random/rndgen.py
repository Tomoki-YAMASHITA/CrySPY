#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import subprocess

from pymatgen import Structure

from . import gen_cell
from . import with_spg
from . import wo_spg


def rndgen_wo_spg(nstruc, natot, atype, nat, cID=0, minlen=4, maxlen=10, dangle=20, mindist=1.5,
                  maxcnt=500, symtoleI=0.001, init_pos_path='./init_POSCARS'):
    '''
    Randomly generate structures without space group information
    '''

    #---------- initialize
    init_struc_data = []
    spgnum = 0

    #---------- cd gen_struc
    if not os.path.isdir('gen_struc'):
        os.mkdir('gen_struc')
    os.chdir('gen_struc')

    #---------- generate structures
    while len(init_struc_data) < nstruc:
        spg_in, a, b, c, alpha, beta, gamma = gen_cell.rndgen_lattice(spgnum, minlen, maxlen, dangle)
        va, vb, vc = gen_cell.calc_latvec(a, b, c, alpha, beta, gamma)
        incoord = wo_spg.gen_coordinate.rndgen_coord(natot, va, vb, vc, mindist, maxcnt)
        if incoord:    # success of generation
            wo_spg.write_poscar.write_poscar(va, vb, vc, incoord, atype, nat)

            #----- pymatgen
            tmp_struc = Structure.from_file('POSCAR')
            if not natot == tmp_struc.num_sites:
                os.remove('POSCAR')
                continue

            #----- check actual space group using pymatgen
            spg_sym, spg_num = tmp_struc.get_space_group_info(symprec=symtoleI)

            #----- register the structure in pymatgen format
            init_struc_data.append(tmp_struc)
            print('Structure ID {0:>8} was generated. Space group: {1:>3} --> {2:>3} {3}'.format(
                   len(init_struc_data) - 1 + cID, spg_in, spg_num, spg_sym))

            #----- save poscar
            save_init_poscar(tmp_struc, len(init_struc_data) - 1 + cID, init_pos_path)

            #----- clear
            os.remove('POSCAR')

    #---------- go back to ..
    os.chdir('../')

    return init_struc_data


def rndgen_spg(nstruc, natot, atype, nat, spgnum='all', cID=0,
               minlen=4, maxlen=10, dangle=20, mindist=1.5,
               maxcnt=500, symtoleI=0.001,
               init_pos_path='./init_POSCARS', fwpath='./find_wy'):
    '''
    Randomly generate structures with space group information
    '''

    #---------- initialize
    init_struc_data = []

    #---------- cd gen_struc
    if not os.path.isdir('gen_struc'):
        os.mkdir('gen_struc')
    os.chdir('gen_struc')

    #---------- generate structures
    while len(init_struc_data) < nstruc:
        spg_in, a, b, c, alpha, beta, gamma = gen_cell.rndgen_lattice(spgnum, minlen, maxlen, dangle)
        cosa, cosb, cosg = gen_cell.calc_cos(alpha, beta, gamma)
        with_spg.fw.fw_input(atype, nat, spg_in, a, b, c, cosa, cosb, cosg)

        #----- loop for same fw_input
        cnt = 0
        while cnt <= maxcnt:
            #----- execute find_wy
            with open('sublog', 'w') as f:
                subprocess.call([fwpath, 'input'], stdout=f, stderr=f)

            #----- generate a structure using POS_WY_SKEL_ALL.json
            if not os.path.isfile('POS_WY_SKEL_ALL.json'):
                wyflag = False
                break
            wyflag, tmp_struc = with_spg.fw.gen_wypos(mindist, maxcnt)
            if wyflag is False:    # Failure
                os.remove('POS_WY_SKEL_ALL.json')
                cnt += 1
                continue
            else:    # Success
                #----- clean
                for rfile in ['input', 'POS_WY_SKEL_ALL.json']:
                    if os.path.isfile(rfile):
                        os.remove(rfile)
                #------ # break fw_input loop
                break

        if wyflag is False:    # maximum trial or no POS_WY_SKEL_ALL.json file
            #----- clean
            for rfile in ['input', 'POS_WY_SKEL_ALL.json']:
                if os.path.isfile(rfile):
                    os.remove(rfile)
            #----- to new fw_input
            continue

        #----- check actual space group using pymatgen
        spg_sym, spg_num = tmp_struc.get_space_group_info(symprec=symtoleI)

        #----- register the structure in pymatgen format
        init_struc_data.append(tmp_struc)
        print('Structure ID {0:>8} was generated. Space group: {1:>3} --> {2:>3} {3}'.format(
               len(init_struc_data) - 1 + cID, spg_in, spg_num, spg_sym))

        #----- save poscar
        save_init_poscar(tmp_struc, len(init_struc_data) - 1 + cID, '../data/init_POSCARS')

        #----- clean
        for rfile in ['input', 'POS_WY_SKEL_ALL.json']:
            if os.path.isfile(rfile):
                os.remove(rfile)

    #---------- go back to ..
    os.chdir('../')

    return init_struc_data


def save_init_poscar(structure, struc_id, fpath):
    #---------- Cart to Direct
    structure.to(fmt='poscar', filename='POSCAR')

    #---------- Change the title of POSCAR
    with open('POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = 'ID_{}\n'.format(struc_id)

    #---------- save POSCAR
    with open(fpath, 'a') as f:
        for line in lines:
            f.write(line)
