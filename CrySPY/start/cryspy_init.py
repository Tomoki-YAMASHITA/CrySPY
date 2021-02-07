'''
Initialize CrySPY
'''

import os

import pandas as pd

from .. import utility
from ..BO import bo_init
from ..EA import ea_init
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..gen_struc.random.gen_pyxtal import Rnd_struc_gen_pyxtal
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

    # ---------- read input
    print('Read input file, cryspy.in')
    rin.readin()          # read input data, cryspy,in
    stat = io_stat.stat_init()    # initialize stat
    rin.writeout()        # write input data in output file, cryspy.out
    rin.save_stat(stat)   # save input variables in cryspy.stat

    # ---------- make data directory
    os.makedirs('data/pkl_data', exist_ok=True)

    # ---------- generate initial structures
    if not rin.load_struc_flag:
        # ------ from scratch
        print('\n# --------- Generate initial structures')
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Generate initial structures\n')
        # ------ pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            rsgx = Rnd_struc_gen_pyxtal(natot=rin.natot, atype=rin.atype,
                                        nat=rin.nat, vol_factor=rin.vol_factor,
                                        vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                                        mindist=rin.mindist,
                                        spgnum=rin.spgnum, symprec=rin.symprec)
            # ------ crystal
            if rin.struc_mode == 'crystal':
                rsgx.gen_struc(nstruc=rin.tot_struc, id_offset=0,
                               init_pos_path='./data/init_POSCARS')
            # ------ molecular crystal
            elif rin.struc_mode == 'mol':
                rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
                rsgx.gen_struc_mol(nstruc=rin.tot_struc, id_offset=0,
                                   init_pos_path='./data/init_POSCARS',
                                   timeout_mol=rin.timeout_mol)
            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
                rsgx.gen_struc_mol_break_sym(nstruc=rin.tot_struc,
                                             rot_mol=rin.rot_mol,
                                             id_offset=0,
                                             init_pos_path='./data/init_POSCARS')
            # ------ init_struc_data
            init_struc_data = rsgx.init_struc_data
        # ------ Rnd_struc_gen
        else:
            rsg = Rnd_struc_gen(natot=rin.natot, atype=rin.atype, nat=rin.nat,
                                minlen=rin.minlen, maxlen=rin.maxlen,
                                dangle=rin.dangle, mindist=rin.mindist,
                                vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                                maxcnt=rin.maxcnt, symprec=rin.symprec)
            if rin.spgnum == 0:
                rsg.gen_wo_spg(nstruc=rin.tot_struc, id_offset=0,
                               init_pos_path='./data/init_POSCARS')
                init_struc_data = rsg.init_struc_data
            else:
                fwpath = utility.check_fwpath()
                rsg.gen_with_find_wy(nstruc=rin.tot_struc, spgnum=rin.spgnum,
                                     id_offset=0, init_pos_path='./data/init_POSCARS',
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
        kpt_data = {}
        pkl_data.save_kpt(kpt_data)
    if rin.energy_step_flag:
        energy_step_data = {}
        pkl_data.save_energy_step(energy_step_data)
    if rin.struc_step_flag:
        struc_step_data = {}
        pkl_data.save_struc_step(struc_step_data)
    if rin.force_step_flag:
        force_step_data = {}
        pkl_data.save_force_step(force_step_data)
    if rin.stress_step_flag:
        stress_step_data = {}
        pkl_data.save_stress_step(stress_step_data)
