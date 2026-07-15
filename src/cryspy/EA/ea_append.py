'''
Append structures by evolutionary algorithm
'''

from logging import getLogger
import os

import pandas as pd

from .ea_child import child_gen
from .natural_selection import natural_selection
from ..IO import change_input, out_results, pkl_data


logger = getLogger('cryspy')


def append_struc(rin, init_struc_data, rng=None):
    # ---------- append structures by EA
    logger.info('# ---------- Append structures by EA')

    # ---------- load data
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()

    # ---------- fitness
    fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}

    # ---------- natural selection
    logger.info('# ------ natural selection')
    ranking, _, _ = natural_selection(
        fitness=fitness,
        struc_data=opt_struc_data,
        elite_struc=None,
        elite_fitness=None,
        n_fittest=rin.n_fittest,
        fit_reverse=rin.fit_reverse,
        emax_ea=rin.emax_ea,
        emin_ea=rin.emin_ea,
        rng=rng,
    )
    logger.info('ranking without duplication:')
    for cid in ranking:
        logger.info(f'Structure ID {cid:>6}, fitness: {fitness[cid]:>10.5f}')

    # ---------- generate children by EA
    logger.info('# ------ Generate children')
    # init_struc_data will be updated and saved in child_gen function
    init_struc_data, parents, operation = child_gen(
        rin=rin,
        ranking=ranking,
        fittest=fitness,
        struc_data=opt_struc_data,
        init_struc_data=init_struc_data,
        nat_data=None,
        rng=rng,
    )

    # ---------- id_queueing
    # id_queueing is treated after this append_struc function

    # ----------  ea_info
    if os.path.isfile('./data/pkl_data/ea_info.pkl'):
        ea_info = pkl_data.load_ea_info()
        ea_origin = pkl_data.load_ea_origin()
    else:
        # ------ initialize
        # -- ea_info
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
        # -- ea_origin
        ea_origin = pd.DataFrame(columns=[
            'Gen',
            'Struc_ID',
            'Operation',
            'Parent',
        ])
        ea_origin.iloc[:, 0:2] = ea_origin.iloc[:, 0:2].astype(int)
    # ------ register and out ea_info
    ea_info.loc[len(ea_info)] = {
        'Gen':         rin.tot_struc,    # Using tot_struc as generation number for append
        'Population':  rin.n_pop,
        'Crossover':   rin.n_crsov,
        'Permutation': rin.n_perm,
        'Strain':      rin.n_strain,
        'Random':      rin.n_rand,
        'Elite':       0,
        'crs_lat':     rin.crs_lat,
        'slct_func':   rin.slct_func,
    }
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    rows = []
    for cid in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        rows.append([
            rin.tot_struc,
            cid,
            operation[cid],
            parents[cid],
        ])
    tmp_origin = pd.DataFrame(rows, columns=ea_origin.columns)
    ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ random part
    rows = []
    for cid in range(rin.tot_struc + rin.n_pop - rin.n_rand,
                     rin.tot_struc + rin.n_pop):
        rows.append([
            rin.tot_struc,
            cid,
            'random',
            None,
        ])
    tmp_origin = pd.DataFrame(rows, columns=ea_origin.columns)
    ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------  out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_data
    pkl_data.save_ea_info(ea_info)
    pkl_data.save_ea_origin(ea_origin)

    # ---------- change variables in cryspy.in
    next_tot = rin.tot_struc + rin.n_pop
    config = change_input.read_config()
    logger.info('# -- Changed cryspy.in')
    # ------ tot_struc
    change_input.change_input(config, 'basic', 'tot_struc', next_tot)
    logger.info(f'Changed the value of tot_struc in cryspy.in from {rin.tot_struc} to {next_tot}')
    rin.tot_struc = next_tot
    # ------ append_struc_ea: True --> False
    change_input.change_input(config, 'option', 'append_struc_ea', False)
    logger.info('Changed the value of append_struc_ea in cryspy.in from True to False')
    # ------ write and save
    change_input.write_config(config)
    pkl_data.save_input(rin)

    # ---------- return
    return init_struc_data
