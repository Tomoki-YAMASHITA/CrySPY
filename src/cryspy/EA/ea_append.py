'''
Append structures by evolutionary algorithm
'''

from logging import getLogger
import os

import pandas as pd

from .ea_child import child_gen
from .gen_struc_EA.select_parents import Select_parents
from ..IO import change_input, io_stat, out_results, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def append_struc(stat, init_struc_data):
    # ---------- append structures by EA
    logger.info('# ---------- Append structures by EA')

    # ---------- load data
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()

    # ---------- fitness
    fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}

    # ---------- instantiate Seclect_parents class
    logger.info('# ------ select parents')
    sp = Select_parents(opt_struc_data, fitness, None, None, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament()
    else:
        sp.set_roulette()

    # ---------- generate offspring by EA
    logger.info('# ------ Generate structures')
    init_struc_data, eagen = child_gen(sp, init_struc_data)

    # ----------  ea_info
    if os.path.isfile('./data/pkl_data/EA_data.pkl'):
        _, _, ea_info, ea_origin = pkl_data.load_ea_data()
    else:
        # ------ initialize
        # -- ea_info
        ea_info = pd.DataFrame(columns=['Gen', 'Population',
                                        'Crossover', 'Permutation', 'Strain',
                                        'Random', 'Elite',
                                        'crs_lat', 'slct_func'])
        ea_info.iloc[:, 0:7] = ea_info.iloc[:, 0:7].astype(int)
        # -- ea_origin
        ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID',
                                          'Operation', 'Parent'])
        ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)
    # ------ register ea_info
    tmp_info = pd.DataFrame([[rin.tot_struc, rin.n_pop, rin.n_crsov,
                              rin.n_perm, rin.n_strain, rin.n_rand, 0,
                              rin.crs_lat, rin.slct_func]],
                            columns=ea_info.columns)
    ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    # ------ out ea_info
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cid in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.DataFrame([[rin.tot_struc, cid, eagen.operation[cid],
                                    eagen.parents[cid]]], columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ random part
    for cid in range(rin.tot_struc + rin.n_pop - rin.n_rand,
                     rin.tot_struc + rin.n_pop):
        tmp_origin = pd.DataFrame([[rin.tot_struc, cid, 'random', None]],
                                  columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------  out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_data
    ea_data = (None, None, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)

    # ---------- change variables in cryspy.in
    config = change_input.config_read()
    logger.info('# -- Changed cryspy.in')
    # ------ tot_struc
    change_input.change_basic(config, 'tot_struc', rin.tot_struc + rin.n_pop)
    logger.info(f'Changed tot_struc in cryspy.in from {rin.tot_struc} to {rin.tot_struc + rin.n_pop}')
    rin.tot_struc = rin.tot_struc + rin.n_pop
    # ------ append_struc_ea: True --> False
    change_input.change_option(config, 'append_struc_ea', False)
    logger.info('Changed append_struc_ea in cryspy.in from True to Flase')
    # ------ write
    change_input.write_config(config)

    # ---------- status
    io_stat.set_input_common(stat, 'basic', 'tot_struc', rin.tot_struc)
    io_stat.set_input_common(stat, 'option', 'append_struc_ea', False)
    io_stat.write_stat(stat)

    # ---------- return
    return init_struc_data
