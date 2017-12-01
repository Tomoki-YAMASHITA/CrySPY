#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser

from .. import utility
from ..gen_struc.random import rndgen
from ..IO import pkl_data
from ..IO import read_input as rin


def restart():
    print('\n\n')
    print(utility.get_date())
    print(utility.get_version())
    print('Restart cryspy.py\n\n')
    stat = ConfigParser.ConfigParser()
    stat.read('cryspy.stat')

    # ---------- read input and check the change
    rin.readin()
    rin.diffinstat(stat)

    return stat


def append_struc(init_struc_data):
    # ---------- append initial structures
    print('\n# ---------- Append structures')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Append structures\n')
    init_pos_path = utility.get_init_pos_path()
    cID = len(init_struc_data)
    nstruc = rin.tot_struc - cID
    if rin.spgnum == 0:
        init_struc_data += rndgen.rndgen_wo_spg(
                               nstruc, rin.natot, rin.atype, rin.nat, cID,
                               rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                               rin.maxcnt, rin.symtoleI, init_pos_path)
    else:
        fwpath = utility.check_fwpath()
        init_struc_data += rndgen.rndgen_spg(
                              nstruc, rin.natot, rin.atype, rin.nat, rin.spgnum, cID,
                              rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                              rin.maxcnt, rin.symtoleI, init_pos_path, fwpath)

    print('')    # for blank line
    with open('cryspy.out', 'a') as fout:
        fout.write('Generated structures up to ID {}\n\n'.format(len(init_struc_data)-1))

    # ---------- save
    pkl_data.save_init_struc(init_struc_data)

    return init_struc_data
