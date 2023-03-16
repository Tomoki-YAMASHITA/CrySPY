'''
Restart CrySPY
'''

import os

from ..IO import io_stat, pkl_data
from ..IO import read_input as rin
from ..util.utility import get_date, get_version, check_fwpath, backup_cryspy
from ..util.struc_util import set_mindist

# ---------- import later
#from ..RS import rs_restart
#from ..BO import bo_restart
#from ..LAQA import laqa_restart
#from ..EA import ea_append
#from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
#from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen


def restart():
    print('\n\n')
    print(get_date())
    print('CrySPY ' + get_version())
    print('Restart cryspy.py\n\n')

    # ---------- read stat
    stat = io_stat.stat_read()

    # ---------- read input and check the change
    rin.readin()
    rin.diffinstat(stat)

    # ---------- load init_struc_data for appending structures
    init_struc_data = pkl_data.load_init_struc()
    if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
        struc_mol_id = pkl_data.load_struc_mol_id()

    # ---------- append structures
    if len(init_struc_data) < rin.tot_struc:
        backup_cryspy()
        prev_nstruc = len(init_struc_data)
        if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
            # ------ mol in EA
            init_struc_data = append_struc(init_struc_data, struc_mol_id)
        else:
            # normal
            init_struc_data = append_struc(init_struc_data, None)
        # ------ RS
        if rin.algo == 'RS':
            from ..RS import rs_restart
            rs_restart.restart(stat, prev_nstruc)
        # ------ BO
        if rin.algo == 'BO':
            from ..BO import bo_restart
            bo_restart.restart(init_struc_data, prev_nstruc)
        # ------ LAQA
        if rin.algo == 'LAQA':
            from ..LAQA import laqa_restart
            laqa_restart.restart(stat, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()
    elif rin.tot_struc < len(init_struc_data):
        raise ValueError('tot_struc < len(init_struc_data)')

    # ---------- append structures by EA (option)
    if rin.append_struc_ea:
        backup_cryspy()
        from ..EA import ea_append
        prev_nstruc = len(init_struc_data)
        init_struc_data = ea_append.append_struc(stat, init_struc_data)
        # ------ RS
        if rin.algo == 'RS':
            from ..RS import rs_restart
            rs_restart.restart(stat, prev_nstruc)
        # ------ BO
        if rin.algo == 'BO':
            from ..BO import bo_restart
            bo_restart.restart(init_struc_data, prev_nstruc)
        # ------ LAQA
        if rin.algo == 'LAQA':
            from ..LAQA import laqa_restart
            laqa_restart.restart(stat, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- return
    return stat, init_struc_data


def append_struc(init_struc_data, struc_mol_id):
    # ---------- append initial structures
    print('\n# ---------- Append structures')
    id_offset = len(init_struc_data)
    nstruc = rin.tot_struc - id_offset

    # ---------- mindist
    print('# ------ mindist')
    mindist = set_mindist(rin.mindist, rin.mindist_factor)

    # ---------- pyxtal
    if not (rin.spgnum == 0 or rin.use_find_wy):
        from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
        rsgx = Rnd_struc_gen_pyxtal(mindist=mindist)
        # ------ crystal
        if rin.struc_mode == 'crystal':
            rsgx.gen_struc(nstruc=nstruc, id_offset=id_offset,
                           init_pos_path='./data/init_POSCARS')
        # ------ molecular crystal
        elif rin.struc_mode == 'mol':
            rsgx.set_mol()
            rsgx.gen_struc_mol(nstruc=nstruc, id_offset=id_offset,
                               init_pos_path='./data/init_POSCARS')
        # ------ molecular crystal breaking symmetry
        elif rin.struc_mode == 'mol_bs':
            print('# -- mindist_mol_bs')
            mindist_dummy = set_mindist(rin.mindist_mol_bs, rin.mindist_mol_bs_factor, dummy=True)
            rsgx.set_mol()
            rsgx.gen_struc_mol_break_sym(nstruc=nstruc,
                                         mindist_dummy=mindist_dummy,
                                         id_offset=id_offset,
                                         init_pos_path='./data/init_POSCARS')
        # ------ update
        init_struc_data.update(rsgx.init_struc_data)
        if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
            struc_mol_id.update(rsgx.struc_mol_id)
    # ---------- w/o pyxtal
    else:
        from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen
        rsg = Rnd_struc_gen(mindist=mindist)
        if rin.spgnum == 0:
            rsg.gen_wo_spg(nstruc=nstruc, id_offset=id_offset,
                           init_pos_path='./data/init_POSCARS')
            init_struc_data.update(rsg.init_struc_data)
        else:
            fwpath = check_fwpath(rin.fwpath)
            rsg.gen_with_find_wy(nstruc=nstruc,
                                 id_offset=id_offset,
                                 init_pos_path='./data/init_POSCARS',
                                 fwpath=fwpath)
            init_struc_data.update(rsg.init_struc_data)

    # ---------- output
    print('')    # for blank line

    # ---------- save
    pkl_data.save_init_struc(init_struc_data)
    if rin.algo == 'EA' and rin.struc_mode in ['mol', 'mol_bs']:
        pkl_data.save_struc_mol_id(struc_mol_id)

    return init_struc_data
