#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser

from .. import utility
from ..gen_struc.random.random_generation import Rnd_struc_gen
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
    id_offset = len(init_struc_data)
    nstruc = rin.tot_struc - id_offset
    rsg = Rnd_struc_gen(rin.natot, rin.atype, rin.nat,
                        rin.minlen, rin.maxlen, rin.dangle,
                        rin.mindist, rin.maxcnt, rin.symprec)
    if rin.spgnum == 0:
        rsg.gen_wo_spg(nstruc, id_offset, init_pos_path='./data/init_POSCARS')
        init_struc_data.update(rsg.init_struc_data)
    else:
        fwpath = utility.check_fwpath()
        rsg.gen_with_spg(nstruc, rin.spgnum, id_offset,
                         init_pos_path='./data/init_POSCARS', fwpath=fwpath)
        init_struc_data.update(rsg.init_struc_data)

    print('')    # for blank line
    with open('cryspy.out', 'a') as fout:
        fout.write('Generated structures up to ID {}\n\n'.format(len(init_struc_data)-1))

    # ---------- save
    pkl_data.save_init_struc(init_struc_data)

    return init_struc_data
