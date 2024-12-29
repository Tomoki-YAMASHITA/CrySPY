from logging import getLogger
import os

import pandas as pd

from ..IO.out_results import out_ea_info, out_nat_data
from ..IO import io_stat, pkl_data
from ..util.struc_util import get_nat


logger = getLogger('cryspy')

def initialize(rin, init_struc_data, rslt_data):
    # ---------- log
    logger.info('# ---------- Initialize evolutionary algorithm')
    logger.info('# ------ Generation 1')
    logger.info(f'{rin.n_pop} structures by random')

    # ---------- initialize
    gen = 1
    id_queueing = [i for i in range(rin.n_pop)]
    id_running = []
    # ------ ea_info
    if rin.algo == 'EA':
        ea_info = pd.DataFrame(columns=[
                                    'Gen',
                                    'Population',
                                    'Crossover',
                                    'Permutation',
                                    'Strain',
                                    'Random',
                                    'Elite',
                                    'crs_lat',
                                    'slct_func',
                                ])
        ea_info.iloc[:, 0:7] = ea_info.iloc[:, 0:7].astype(int)
        tmp_info = pd.DataFrame([[
                                    1,
                                    rin.n_pop,
                                    0,
                                    0,
                                    0,
                                    rin.n_pop,
                                    0,
                                    rin.crs_lat,
                                    rin.slct_func,
                                ]], columns=ea_info.columns)
        ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    if rin.algo == 'EA-vc':
        ea_info = pd.DataFrame(columns=[
                                    'Gen',
                                    'Population',
                                    'Crossover',
                                    'Permutation',
                                    'Strain',
                                    'Addition',
                                    'Elimination',
                                    'Substitution',
                                    'Random',
                                    'Elite',
                                    'crs_lat',
                                    'slct_func',
                                ])
        ea_info.iloc[:, 0:10] = ea_info.iloc[:, 0:10].astype(int)
        tmp_info = pd.DataFrame([[
                            1,
                            rin.n_pop,
                            0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            rin.n_pop,
                            0,
                            rin.crs_lat,
                            rin.slct_func,
                    ]], columns=ea_info.columns)
        ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    out_ea_info(ea_info)
    # ------ ea_origin
    ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID', 'Operation', 'Parent'])
    ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)
    for cid in range(rin.n_pop):
        tmp_origin = pd.DataFrame([[1, cid, 'random', None]],
                                  columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ elite
    elite_struc = None
    elite_fitness = None
    # ------ rslt_data
    rslt_data['Gen'] = pd.Series(dtype=int)
    rslt_data = rslt_data[[
                    'Gen',
                    'Spg_num',
                    'Spg_sym',
                    'Spg_num_opt',
                    'Spg_sym_opt',
                    'E_eV_atom',
                    'Magmom',
                    'Opt',
                ]]
    if rin.algo == 'EA-vc':
        rslt_data['Ef_eV_atom'] = pd.Series(dtype='float64')
        rslt_data['Num_atom'] = pd.Series(dtype='object')
        rslt_data = rslt_data[[
                        'Gen',
                        'Spg_num',
                        'Spg_sym',
                        'Spg_num_opt',
                        'Spg_sym_opt',
                        'E_eV_atom',
                        'Ef_eV_atom',
                        'Num_atom',
                        'Magmom',
                        'Opt',
                    ]]
    # ------ nat, hdist, etc for EA-vc
    if rin.algo == 'EA-vc':
        nat_data = {}      # {ID: [nat], ...}
        hdist_data = {}    # {gen: {ID: hdist, ...}, ...}
        pd_data = {}       # {gen: pd, ...}
        for cid, struc in init_struc_data.items():
            tmp_nat = get_nat(struc, rin.atype)
            nat_data[cid] = tmp_nat
        out_nat_data(nat_data, rin.atype)
        os.makedirs('data/convex_hull', exist_ok=True)

    # ---------- save
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)
    pkl_data.save_gen(gen)
    pkl_data.save_elite_struc(elite_struc)
    pkl_data.save_elite_fitness(elite_fitness)
    pkl_data.save_ea_info(ea_info)
    pkl_data.save_ea_origin(ea_origin)
    pkl_data.save_rslt(rslt_data)
    if rin.algo == 'EA-vc':
        pkl_data.save_nat_data(nat_data)
        pkl_data.save_hdist_data(hdist_data)
        pkl_data.save_pd_data(pd_data)

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_common(stat, 'generation', gen)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)
