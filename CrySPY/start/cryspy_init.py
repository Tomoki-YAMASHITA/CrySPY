'''
Initialize CrySPY
'''

import os

import pandas as pd

from .. import utility
from ..BO import bo_init
from ..EA import ea_init
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..IO import pkl_data, io_stat
from ..IO import read_input as rin
from ..LAQA import laqa_init
from ..RS import rs_init


def initialize():
    # ---------- start
    print(utility.get_date())
    print(utility.get_version())
    print('Start cryspy.py\n')
    with open('cryspy.out', 'w') as fout:
        fout.write(utility.get_date() + '\n')
        fout.write(utility.get_version() + '\n')
        fout.write('Start cryspy.py\n\n')

    # ---------- initialize stat
    stat = io_stat.stat_init()

    # ---------- read input
    print('Read input file, cryspy.in')
    rin.readin()          # read input data, cryspy,in
    rin.writeout()        # write input data in output file, cryspy.out
    rin.save_stat(stat)   # save input variables in cryspy.stat

    # ---------- make data directory
    if not os.path.isdir('data/pkl_data'):
        print('Make data directory ./data/pkl_data')
        os.makedirs('data/pkl_data')

    # ---------- generate initial structures
    if not rin.load_struc_flag:
        # ------ from scratch
        print('\n# --------- Generate initial structures')
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Generate initial structures\n')
        rsg = Rnd_struc_gen(rin.natot, rin.atype, rin.nat,
                            rin.minlen, rin.maxlen, rin.dangle,
                            rin.mindist, rin.maxcnt, rin.symprec)
        if rin.spgnum == 0:
            rsg.gen_wo_spg(rin.tot_struc, id_offset=0,
                           init_pos_path='./data/init_POSCARS')
            init_struc_data = rsg.init_struc_data
        else:
            fwpath = utility.check_fwpath()
            rsg.gen_with_spg(rin.tot_struc, rin.spgnum, id_offset=0,
                             init_pos_path='./data/init_POSCARS',
                             fwpath=fwpath)
            init_struc_data = rsg.init_struc_data
        with open('cryspy.out', 'a') as fout:
            fout.write('Generated structures up to ID {}\n\n'.format(
                len(init_struc_data)-1))
        # ------ save
        pkl_data.save_init_struc(init_struc_data)
    else:
        # ------ load initial structure
        print('\n# --------- Load initial structure data')
        print('Load ./data/pkl_data/init_struc_data.pkl\n')
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Load initial structure data\n')
            fout.write('Load ./data/pkl_data/init_struc_data.pkl\n\n')
        init_struc_data = pkl_data.load_init_struc()
        # -- check
        if not rin.tot_struc == len(init_struc_data):
            raise ValueError('rin.tot_struc = {0},'
                             ' len(init_struc_data) = {1}'.format(
                                 rin.tot_struc, len(init_struc_data)))

    # ---------- initialize opt_struc_data
    opt_struc_data = {}
    pkl_data.save_opt_struc(opt_struc_data)

    # ---------- initialize rslt_data
    rslt_data = pd.DataFrame(columns=['Spg_num', 'Spg_sym',
                                      'Spg_num_opt', 'Spg_sym_opt',
                                      'E_eV_atom', 'Magmom', 'Opt'])
    rslt_data[['Spg_num', 'Spg_num_opt']] = rslt_data[
                                   ['Spg_num', 'Spg_num_opt']].astype(int)
    pkl_data.save_rslt(rslt_data)

    # ---------- initialize for each algorithm
    if rin.algo == 'RS':
        rs_init.initialize(stat)
    elif rin.algo == 'BO':
        bo_init.initialize(stat, init_struc_data, rslt_data)
    elif rin.algo == 'LAQA':
        laqa_init.initialize(stat)
    elif rin.algo == "EA":
        ea_init.initialize(stat, rslt_data)

    # ---------- initialize etc
    if rin.kpt_flag:
        kpt_init()
    if rin.energy_step_flag:
        energy_step_init()
    if rin.struc_step_flag:
        struc_step_init()
    if rin.fs_step_flag:
        fs_step_init()


def kpt_init():
    kpt_data = {}
    pkl_data.save_kpt(kpt_data)


def energy_step_init():
    energy_step_data = {}
    pkl_data.save_energy_step(energy_step_data)


def struc_step_init():
    struc_step_data = {}
    pkl_data.save_struc_step(struc_step_data)


def fs_step_init():
    force_step_data = {}
    stress_step_data = {}
    fs_step_data = (force_step_data, stress_step_data)
    pkl_data.save_fs_step(fs_step_data)
