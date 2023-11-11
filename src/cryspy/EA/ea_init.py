'''
Initialize evolutionary algorithm
'''

from logging import getLogger

import pandas as pd

from ..IO import out_results
from ..IO import io_stat, pkl_data
from ..IO import read_input as rin
from ..util.struc_util import get_nat


logger = getLogger('cryspy')

def initialize(stat, init_struc_data, rslt_data):
    # ---------- log
    logger.info('# ---------- Initialize evolutionary algorithm')
    logger.info('# ------ Generation 1')
    logger.info(f'{rin.tot_struc} structures by random')

    # ---------- initialize
    gen = 1
    id_queueing = [i for i in range(rin.tot_struc)]
    id_running = []
    # ------ ea_info
    ea_info = pd.DataFrame(columns=['Gen', 'Population',
                                    'Crossover', 'Permutation', 'Strain',
                                    'Random', 'Elite',
                                    'crs_lat', 'slct_func'])
    ea_info.iloc[:, 0:7] = ea_info.iloc[:, 0:7].astype(int)
    tmp_info = pd.DataFrame([[1, rin.tot_struc, 0, 0, 0, rin.tot_struc, 0,
                              rin.crs_lat, rin.slct_func]],
                            columns=ea_info.columns)
    ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    out_results.out_ea_info(ea_info)
    # ------ ea_origin
    ea_origin = pd.DataFrame(columns=['Gen', 'Struc_ID',
                                      'Operation', 'Parent'])
    ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)
    for cid in range(rin.tot_struc):
        tmp_origin = pd.DataFrame([[1, cid, 'random', None]],
                                  columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ elite
    elite_struc = None
    elite_fitness = None
    # ------ rslt_data
    rslt_data['Gen'] = pd.Series(dtype=int)
    rslt_data = rslt_data[['Gen', 'Spg_num',
                           'Spg_sym', 'Spg_num_opt',
                           'Spg_sym_opt', 'E_eV_atom', 'Magmom', 'Opt']]
    if rin.algo == 'EA-vc':
        rslt_data['Ef_eV_atom'] = pd.Series(dtype='float64')
        rslt_data = rslt_data[['Gen', 'Spg_num', 'Spg_sym',
                               'Spg_num_opt', 'Spg_sym_opt',
                               'E_eV_atom', 'Ef_eV_atom', 'Magmom', 'Opt']]
    # ------ nat, ratio for EA-vc
    if rin.algo == 'EA-vc':
        nat_data = {}
        ratio_data = {}
        for cid, struc in init_struc_data.items():
            tmp_nat, tmp_ratio = get_nat(struc, rin.atype)
            nat_data[cid] = tmp_nat
            ratio_data[cid] = tmp_ratio

    # ---------- save
    ea_id_data = (gen, id_queueing, id_running)
    pkl_data.save_ea_id(ea_id_data)
    ea_data = (elite_struc, elite_fitness, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)
    pkl_data.save_rslt(rslt_data)
    if rin.algo == 'EA-vc':
        ea_vc_data = (nat_data, ratio_data)
        pkl_data.save_ea_vc_data(ea_vc_data)

    # ---------- status
    io_stat.set_common(stat, 'generation', gen)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)
