'''
Restart CrySPY
'''

import os

from .. import utility
from ..BO import bo_restart
from ..EA import ea_append
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..gen_struc.random.gen_pyxtal import Rnd_struc_gen_pyxtal
from ..IO import io_stat, pkl_data
from ..IO import read_input as rin
from ..LAQA import laqa_restart
from ..RS import rs_restart


def restart():
    print('\n\n')
    print(utility.get_date())
    print(utility.get_version())
    print('Restart cryspy.py\n\n')

    # ---------- read stat
    stat = io_stat.stat_read()

    # ---------- read input and check the change
    rin.readin()
    rin.diffinstat(stat)

    # ---------- load init_struc_data for appending structures
    init_struc_data = pkl_data.load_init_struc()

    # ---------- append structures
    if len(init_struc_data) < rin.tot_struc:
        prev_nstruc = len(init_struc_data)
        init_struc_data = append_struc(init_struc_data)
        # ------ RS
        if rin.algo == 'RS':
            rs_restart.restart(stat, prev_nstruc)
        # ------ BO
        if rin.algo == 'BO':
            bo_restart.restart(init_struc_data, prev_nstruc)
        # ------ LAQA
        if rin.algo == 'LAQA':
            laqa_restart.restart(stat, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()
    elif rin.tot_struc < len(init_struc_data):
        raise ValueError('tot_struc < len(init_struc_data)')

    # ---------- append structures by EA (option)
    if rin.append_struc_ea:
        prev_nstruc = len(init_struc_data)
        init_struc_data = ea_append.append_struc(stat, init_struc_data)
        # ------ RS
        if rin.algo == 'RS':
            rs_restart.restart(stat, prev_nstruc)
        # ------ BO
        if rin.algo == 'BO':
            bo_restart.restart(init_struc_data, prev_nstruc)
        # ------ LAQA
        if rin.algo == 'LAQA':
            laqa_restart.restart(stat, prev_nstruc)
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- return
    return stat, init_struc_data


def append_struc(init_struc_data):
    # ---------- append initial structures
    print('\n# ---------- Append structures')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Append structures\n')
    id_offset = len(init_struc_data)
    nstruc = rin.tot_struc - id_offset

    # ---------- pyxtal
    if not (rin.spgnum == 0 or rin.use_find_wy):
        rsgx = Rnd_struc_gen_pyxtal(natot=rin.natot, atype=rin.atype,
                                    nat=rin.nat, vol_factor=rin.vol_factor,
                                    vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                                    mindist=rin.mindist,
                                    spgnum=rin.spgnum, symprec=rin.symprec)
        # ------ crystal
        if rin.struc_mode == 'crystal':
            rsgx.gen_struc(nstruc=nstruc, id_offset=id_offset,
                           init_pos_path='./data/init_POSCARS')
        # ------ molecular crystal
        elif rin.struc_mode == 'mol':
            rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
            rsgx.gen_struc_mol(nstruc=nstruc, id_offset=id_offset,
                               init_pos_path='./data/init_POSCARS',
                               timeout_mol=rin.timeout_mol)
        # ------ molecular crystal breaking symmetry
        elif rin.struc_mode == 'mol_bs':
            rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
            rsgx.gen_struc_mol_break_sym(nstruc=nstruc,
                                         rot_mol=rin.rot_mol,
                                         id_offset=id_offset,
                                         init_pos_path='./data/init_POSCARS')
        # ------ update
        init_struc_data.update(rsgx.init_struc_data)
    # ---------- Rnd_struc_gen
    else:
        rsg = Rnd_struc_gen(natot=rin.natot, atype=rin.atype, nat=rin.nat,
                            minlen=rin.minlen, maxlen=rin.maxlen,
                            dangle=rin.dangle, mindist=rin.mindist,
                            vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                            maxcnt=rin.maxcnt, symprec=rin.symprec)
        if rin.spgnum == 0:
            rsg.gen_wo_spg(nstruc=nstruc, id_offset=id_offset,
                           init_pos_path='./data/init_POSCARS')
            init_struc_data.update(rsg.init_struc_data)
        else:
            fwpath = utility.check_fwpath()
            rsg.gen_with_find_wy(nstruc=nstruc, spgnum=rin.spgnum,
                                 id_offset=id_offset,
                                 init_pos_path='./data/init_POSCARS',
                                 fwpath=fwpath)
            init_struc_data.update(rsg.init_struc_data)

    # ---------- output
    print('')    # for blank line
    with open('cryspy.out', 'a') as fout:
        fout.write('Generated structures up to ID {}\n\n'.format(
            len(init_struc_data)-1))

    # ---------- save
    pkl_data.save_init_struc(init_struc_data)

    return init_struc_data
