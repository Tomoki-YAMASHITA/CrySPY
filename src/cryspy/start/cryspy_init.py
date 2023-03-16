'''
Initialize CrySPY
'''

import os

import pandas as pd

from ..IO import pkl_data, io_stat
from ..IO import read_input as rin
from ..util.utility import get_date, get_version, check_fwpath
from ..util.struc_util import set_mindist

# ---------- import later
#from ..RS import rs_init
#from ..BO import bo_init
#from ..LAQA import laqa_init
#from ..EA import ea_init
#from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
#from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen


def initialize():
    # ---------- start
    print(get_date())
    print('CrySPY ' + get_version())
    print('Start cryspy.py\n', flush=True)

    # ---------- read input
    print('Read input file, cryspy.in')
    rin.readin()          # read input data, cryspy,in
    stat = io_stat.stat_init()    # initialize stat
    rin.save_stat(stat)   # save input variables in cryspy.stat

    # ---------- make data directory
    os.makedirs('data/pkl_data', exist_ok=True)

    # ---------- generate initial structures
    if not rin.load_struc_flag:
        # ------ from scratch
        print('\n# --------- Generate initial structures')
        # ------ mindist
        print('# ------ mindist')
        mindist = set_mindist(rin.mindist, rin.mindist_factor)
        # ------ pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
            rsgx = Rnd_struc_gen_pyxtal(mindist=mindist)
            # ------ crystal
            if rin.struc_mode == 'crystal':
                rsgx.gen_struc(nstruc=rin.tot_struc, id_offset=0,
                               init_pos_path='./data/init_POSCARS')
            # ------ molecular crystal
            elif rin.struc_mode == 'mol':
                rsgx.set_mol()
                rsgx.gen_struc_mol(nstruc=rin.tot_struc, id_offset=0,
                                   init_pos_path='./data/init_POSCARS')
            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                print('# -- mindist_mol_bs')
                mindist_dummy = set_mindist(rin.mindist_mol_bs, rin.mindist_mol_bs_factor, dummy=True)
                rsgx.set_mol()
                rsgx.gen_struc_mol_break_sym(nstruc=rin.tot_struc,
                                             mindist_dummy=mindist_dummy,
                                             id_offset=0,
                                             init_pos_path='./data/init_POSCARS')
            # ------ init_struc_data
            init_struc_data = rsgx.init_struc_data
            if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
                struc_mol_id = rsgx.struc_mol_id
        # ------ w/o pyxtal
        else:
            from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen
            rsg = Rnd_struc_gen(mindist=mindist)
            if rin.spgnum == 0:
                rsg.gen_wo_spg(nstruc=rin.tot_struc, id_offset=0,
                               init_pos_path='./data/init_POSCARS')
                init_struc_data = rsg.init_struc_data
            else:
                fwpath = check_fwpath(rin.fwpath)
                rsg.gen_with_find_wy(nstruc=rin.tot_struc,
                                     id_offset=0, init_pos_path='./data/init_POSCARS',
                                     fwpath=fwpath)
                init_struc_data = rsg.init_struc_data
        # ------ save
        pkl_data.save_init_struc(init_struc_data)
        if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
            pkl_data.save_struc_mol_id(struc_mol_id)
    else:
        # ------ load initial structure
        print('\n# --------- Load initial structure data')
        print('Load ./data/pkl_data/init_struc_data.pkl\n')
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
        from ..RS import rs_init
        rs_init.initialize(stat)
    elif rin.algo == 'BO':
        from ..BO import bo_init
        bo_init.initialize(stat, init_struc_data, rslt_data)
    elif rin.algo == 'LAQA':
        from ..LAQA import laqa_init
        laqa_init.initialize(stat)
    elif rin.algo == "EA":
        from ..EA import ea_init
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

    # ---------- for ext
    if rin.calc_code == 'ext':
        os.makedirs('ext', exist_ok=True)
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('out\n')
