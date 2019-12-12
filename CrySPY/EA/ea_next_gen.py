#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser
import pandas as pd

from .. import utility
from ..gen_struc.EA.select_parents import Select_parents
from ..gen_struc.EA.crossover import Crossover
from ..gen_struc.EA.permutation import Permutation
from ..gen_struc.EA.strain import Strain
from ..gen_struc.EA.ea_generation import EA_generation
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..IO import out_results
from ..IO import pkl_data
from ..IO import read_input as rin


def next_gen(stat, init_struc_data, opt_struc_data, rslt_data, ea_id_data):
    # ---------- ea_id_data
    gen, next_id = ea_id_data

    # ---------- out and log
    with open('cryspy.out', 'a') as fout:
        fout.write('# ---------- Evolutionary algorithm\n\n')
        fout.write('# ------ Generation {}\n'.format(gen + 1))
    print('# ---------- Evolutionary algorithm\n')
    print('# ------ Generation {}'.format(gen + 1))

    # ---------- current generation
    c_rslt = rslt_data[rslt_data['Gen'] == gen]
    c_fitness = dict(zip(c_rslt['Struc_ID'].values, c_rslt['E_eV_atom'].values))

    # ---------- load ea_data, ea_data is used only in this module
    elite_struc, elite_fitness, ea_info, ea_origin = pkl_data.load_ea_data()

    # ---------- instantiate Seclect_parents class
    print('# -- select parents')
    sp = Select_parents(opt_struc_data, c_fitness, elite_struc, elite_fitness,
                        rin.fit_reverse, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament(t_size=rin.t_size)
    else:
        sp.set_roulette(a=rin.a_rlt, b=rin.b_rlt)

    # ---------- generate offspring by EA
    print('# -- Generate structures')
    eagen = EA_generation(sp=sp, symprec=rin.symprec, id_start=rin.tot_struc,
                          init_pos_path='./data/init_POSCARS')
    # ------ instantiate Crossover class
    if rin.n_crsov > 0:
        co = Crossover(rin.atype, rin.nat, rin.mindist, rin.crs_lat, rin.crs_func,
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
            rsg.gen_wo_spg(rin.n_rand, id_offset=eagen.cID, init_pos_path='./data/init_POSCARS')
            init_struc_data.update(rsg.init_struc_data)
        else:
            fwpath = utility.check_fwpath()
            rsg.gen_with_spg(rin.n_rand, rin.spgnum, id_offset=eagen.cID,
                             init_pos_path='./data/init_POSCARS', fwpath=fwpath)
            init_struc_data.update(rsg.init_struc_data)
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by random\n'.format(rin.n_rand))

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)

    # ---------- select elite
    if rin.n_elite > 0:
        print('# -- select elites')
        # ------ init
        all_fitness = dict(zip(rslt_data['Struc_ID'].values, rslt_data['E_eV_atom'].values))
        elite_struc = {}
        elite_fitness = {}
        # ------ Select_parents class also works for selecting elite structures
        se = Select_parents(opt_struc_data, all_fitness, None, None,
                            rin.fit_reverse, rin.n_elite)
        for cID in se.ranking_dedupe:
            print('Structure ID {0:>8} keeps as the elite'.format(cID))
            elite_struc[cID] = opt_struc_data[cID]
            elite_fitness[cID] = all_fitness[cID]
    else:
        elite_struc = None
        elite_fitness = None
    # ---------- out
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures keeps as the elite\n'.format(rin.n_elite))

    # ---------- new gen
    gen += 1

    # ---------- ea_info
    tmp_info = pd.Series([gen, rin.n_pop, rin.n_crsov, rin.n_perm,
                          rin.n_strain, rin.n_rand, rin.n_elite,
                          rin.crs_func, rin.crs_lat, rin.slct_func], index=ea_info.columns)
    ea_info = ea_info.append(tmp_info, ignore_index=True)
    # ------ out ea_info
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cID in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.Series([gen, cID, eagen.operation[cID],
                                eagen.parents[cID]], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ random part
    for cID in range(rin.tot_struc + rin.n_pop - rin.n_rand, rin.tot_struc + rin.n_pop):
        tmp_origin = pd.Series([gen, cID, 'random', None], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ elite part
    if rin.n_elite > 0:
        for cID in se.ranking_dedupe:
            tmp_origin = pd.Series([gen, cID, 'elite', 'elite'], index=ea_origin.columns)
            ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_id_data
    ea_id_data = (gen, next_id)
    pkl_data.save_ea_id(ea_id_data)

    # ---------- save ea_data
    ea_data = (elite_struc, elite_fitness, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)

    # ---------- change the value of tot_struc
    config = ConfigParser.ConfigParser()
    config.read('cryspy.in')
    config.set('basic', 'tot_struc', '{}'.format(rin.tot_struc + rin.n_pop))
    with open('cryspy.in', 'w') as f:
        config.write(f)
    print('# -- changed cryspy.in')
    print('Changed the value of tot_struc in cryspy.in from {} to {}'.format(
          rin.tot_struc, rin.tot_struc + rin.n_pop))

    # ---------- status
    stat.set('input', 'tot_struc', '{}'.format(rin.tot_struc + rin.n_pop))
    stat.set('status', 'generation', '{}'.format(gen))
    stat.set('status', 'next_id', '{}'.format(next_id))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)
