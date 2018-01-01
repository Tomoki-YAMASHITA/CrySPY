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
    id_offset = len(init_struc_data)
    nstruc = rin.tot_struc - id_offset
    if rin.spgnum == 0:
        tmp_struc_data = rndgen.rndgen_wo_spg(
                               nstruc, rin.natot, rin.atype, rin.nat, id_offset,
                               rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                               rin.maxcnt, rin.symtoleI, init_pos_path)
        init_struc_data.update(tmp_struc_data)
    else:
        fwpath = utility.check_fwpath()
        tmp_struc_data = rndgen.rndgen_spg(
                              nstruc, rin.natot, rin.atype, rin.nat, rin.spgnum, id_offset,
                              rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                              rin.maxcnt, rin.symtoleI, init_pos_path, fwpath)
        init_struc_data.update(tmp_struc_data)

    print('')    # for blank line
    with open('cryspy.out', 'a') as fout:
        fout.write('Generated structures up to ID {}\n\n'.format(len(init_struc_data)-1))

    # ---------- save
    pkl_data.save_init_struc(init_struc_data)

    return init_struc_data
