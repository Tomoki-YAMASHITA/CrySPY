#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import ConfigParser
import os

import pandas as pd

from .. import utility
from ..gen_struc.EA.select_parents import Select_parents
from ..gen_struc.EA.crossover import Crossover
from ..gen_struc.EA.permutation import Permutation
from ..gen_struc.EA.strain import Strain
from ..gen_struc.EA.ea_generation import EA_generation
from ..gen_struc.random import rndgen
from ..IO import out_results
from ..IO import pkl_data
from ..IO import read_input as rin


def append_struc(stat, init_struc_data, opt_struc_data, rslt_data):
    # ---------- append structures by EA
    print('\n# ---------- Append structures by EA')
    with open('cryspy.out', 'a') as fout:
        fout.write('\n# ---------- Append structures by EA\n')

    # ---------- fitness
    fitness = dict(zip(rslt_data['Struc_ID'].values, rslt_data['Energy'].values))

    # ---------- instantiate Seclect_parents class
    print('# ------ select parents')
    sp = Select_parents(opt_struc_data, fitness, None, None, rin.fit_reverse, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament(t_size=rin.t_size)
    else:
        sp.set_roulette(a=rin.a_rlt, b=rin.b_rlt)

    # ---------- generate offspring by EA
    print('# ------ Generate structures')
    eagen = EA_generation(sp=sp, symprec=rin.symprec, id_start=rin.tot_struc, init_pos_path='./data/init_POSCARS')
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
        if rin.spgnum == 0:
            tmp_struc_data = rndgen.rndgen_wo_spg(
                                   rin.n_rand, rin.natot, rin.atype, rin.nat, eagen.cID,
                                   rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                                   rin.maxcnt, rin.symprec, '../data/init_POSCARS')
            # ------ update init_struc_data
            init_struc_data.update(tmp_struc_data)
        else:
            fwpath = utility.check_fwpath()
            tmp_struc_data = rndgen.rndgen_spg(
                                  rin.n_rand, rin.natot, rin.atype, rin.nat, rin.spgnum, eagen.cID,
                                  rin.minlen, rin.maxlen, rin.dangle, rin.mindist,
                                  rin.maxcnt, rin.symprec, '../data/init_POSCARS', fwpath)
            # ------ update init_struc_data
            init_struc_data.update(tmp_struc_data)
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
        ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID', 'Operation', 'Parent'])
        ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)

    # ---------- ea_info
    tmp_info = pd.Series([rin.tot_struc, rin.n_pop, rin.n_crsov, rin.n_perm, rin.n_strain, rin.n_rand, 0,
                          rin.crs_func, rin.crs_lat, rin.slct_func], index=ea_info.columns)
    ea_info = ea_info.append(tmp_info, ignore_index=True)
    # ------ out ea_info
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cID in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.Series([rin.tot_struc, cID, eagen.operation[cID], eagen.parents[cID]], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    # ------ random part
    for cID in range(rin.tot_struc + rin.n_pop - rin.n_rand, rin.tot_struc + rin.n_pop):
        tmp_origin = pd.Series([rin.tot_struc, cID, 'random', None], index=ea_origin.columns)
        ea_origin = ea_origin.append(tmp_origin, ignore_index=True)
    #------  out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_data
    ea_data = (None, None, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)

    # ---------- change variables in cryspy.in
    config = ConfigParser.ConfigParser()
    config.read('cryspy.in')
    print('# -- changed cryspy.in')
    # ------ tot_struc
    config.set('basic', 'tot_struc', '{}'.format(rin.tot_struc + rin.n_pop))
    print('Changed the value of tot_struc in cryspy.in from {} to {}'.format(
          rin.tot_struc, rin.tot_struc + rin.n_pop))
    # ------ append_struc_ea
    config.set('option', 'append_struc_ea', '{}'.format(False))
    print('Changed the value of append_struc_ea in cryspy.in from {} to {}'.format(
          True, False))
    # ------ write
    with open('cryspy.in', 'w') as f:
        config.write(f)

    # ---------- status
    stat.set('input', 'tot_struc', '{}'.format(rin.tot_struc + rin.n_pop))
    stat.set('input', 'append_struc_ea', '{}'.format(False))
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)

    # ---------- return
    return init_struc_data
