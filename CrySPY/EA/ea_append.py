'''
Append structures by evolutionary algorithm
'''

import os

import pandas as pd

from .. import utility
from ..gen_struc.EA.select_parents import Select_parents
from ..gen_struc.EA.crossover import Crossover
from ..gen_struc.EA.permutation import Permutation
from ..gen_struc.EA.strain import Strain
from ..gen_struc.EA.ea_generation import EA_generation
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..IO import out_results
from ..IO import change_input, io_stat, pkl_data
from ..IO import read_input as rin


def append_struc(stat, init_struc_data):
    # ---------- append structures by EA
    print('\n# ---------- Append structures by EA')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Append structures by EA\n')

    # ---------- load data
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()

    # ---------- fitness
    fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}

    # ---------- instantiate Seclect_parents class
    print('# ------ select parents')
    sp = Select_parents(opt_struc_data, fitness, None, None,
                        rin.fit_reverse, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament(t_size=rin.t_size)
    else:
        sp.set_roulette(a=rin.a_rlt, b=rin.b_rlt)

    # ---------- generate offspring by EA
    print('# ------ Generate structures')
    eagen = EA_generation(sp=sp, symprec=rin.symprec, id_start=rin.tot_struc,
                          init_pos_path='./data/init_POSCARS')
    # ------ instantiate Crossover class
    if rin.n_crsov > 0:
        co = Crossover(rin.atype, rin.nat, rin.mindist,
                       rin.crs_lat, rin.crs_func,
                       rin.nat_diff_tole, rin.maxcnt_ea)
        eagen.gen_crossover(rin.n_crsov, co=co)    # crossover
        with open('cryspy.out', 'a') as fout:
            fout.write('{} structures by crossover\n'.format(rin.n_crsov))
    # ------ instantiate Permutation class
    if rin.n_perm > 0:
        pm = Permutation(rin.atype, rin.mindist, rin.ntimes, rin.maxcnt_ea)
        eagen.gen_permutation(rin.n_perm, pm=pm)    # permutation
        with open('cryspy.out', 'a') as fout:
            fout.write('{} structures by permutation\n'.format(rin.n_perm))
    # ------ instantiate Strain class
    if rin.n_strain > 0:
        st = Strain(rin.atype, rin.mindist, rin.sigma_st, rin.maxcnt_ea)
        eagen.gen_strain(rin.n_strain, st=st)    # strain
        with open('cryspy.out', 'a') as fout:
            fout.write('{} structures by strain\n'.format(rin.n_strain))
    # ------ update init_struc_data
    init_struc_data.update(eagen.offspring)

    # ---------- random generation
    if rin.n_rand > 0:
        rsg = Rnd_struc_gen(rin.natot, rin.atype, rin.nat,
                            rin.minlen, rin.maxlen, rin.dangle,
                            rin.mindist, rin.maxcnt, rin.symprec)
        if rin.spgnum == 0:
            rsg.gen_wo_spg(rin.n_rand, id_offset=eagen.cid,
                           init_pos_path='./data/init_POSCARS')
            init_struc_data.update(rsg.init_struc_data)
        else:
            fwpath = utility.check_fwpath()
            rsg.gen_with_spg(rin.n_rand, rin.spgnum, id_offset=eagen.cid,
                             init_pos_path='./data/init_POSCARS',
                             fwpath=fwpath)
            init_struc_data.update(rsg.init_struc_data)
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by random\n'.format(rin.n_rand))

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)

    # ---------- load or init ea_data
    if os.path.isfile('./data/pkl_data/EA_data.pkl'):
        _, _, ea_info, ea_origin = pkl_data.load_ea_data()
    else:
        # ------ initialize
        # -- ea_info
        ea_info = pd.DataFrame(columns=['Gen', 'Population',
                                        'Crossover', 'Permutation', 'Strain',
                                        'Random', 'Elite',
                                        'crs_func', 'crs_lat', 'slct_func'])
        ea_info.iloc[:, 0:7] = ea_info.iloc[:, 0:7].astype(int)
        # -- ea_origin
        ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID',
                                          'Operation', 'Parent'])
        ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)

    # ---------- ea_info
    tmp_info = pd.Series([rin.tot_struc, rin.n_pop, rin.n_crsov,
                          rin.n_perm, rin.n_strain, rin.n_rand, 0,
                          rin.crs_func, rin.crs_lat, rin.slct_func],
                         index=ea_info.columns)
    ea_info = ea_info.append(tmp_info, ignore_index=True)
    # ------ out ea_info
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cid in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.Series([rin.tot_struc, cid, eagen.operation[cid],
                                eagen.parents[cid]], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ random part
    for cid in range(rin.tot_struc + rin.n_pop - rin.n_rand,
                     rin.tot_struc + rin.n_pop):
        tmp_origin = pd.Series([rin.tot_struc, cid, 'random', None],
                               index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------  out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_data
    ea_data = (None, None, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)

    # ---------- change variables in cryspy.in
    config = change_input.config_read()
    print('# -- Changed cryspy.in')
    # ------ tot_struc
    change_input.change_basic(config, 'tot_struc', rin.tot_struc + rin.n_pop)
    print('Changed tot_struc in cryspy.in from {} to {}'.format(
          rin.tot_struc, rin.tot_struc + rin.n_pop))
    rin.tot_struc = rin.tot_struc + rin.n_pop
    # ------ append_struc_ea: True --> False
    change_input.change_option(config, 'append_struc_ea', False)
    print('Changed append_struc_ea in cryspy.in from {} to {}'.format(
          True, False))
    # ------ write
    change_input.write_config(config)

    # ---------- status
    io_stat.set_input_common(stat, 'tot_struc', rin.tot_struc + rin.n_pop)
    io_stat.set_input_common(stat, 'append_struc_ea', False)
    io_stat.write_stat(stat)

    # ---------- return
    return init_struc_data
