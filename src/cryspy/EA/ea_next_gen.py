'''
Generational change in evolutionary algorithm
'''

from logging import getLogger

import pandas as pd

from .calc_hull import write_asc_hdist
from .ea_child import child_gen
from .gen_struc_EA.select_parents import Select_parents
from ..IO import out_results
from ..IO import change_input, io_stat, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def next_gen(stat, init_struc_data, struc_mol_id, opt_struc_data, rslt_data, ea_id_data):
    # ---------- ea_id_data
    gen, id_queueing, id_running = ea_id_data

    # ---------- log
    logger.info('# ---------- Evolutionary algorithm')
    logger.info(f'# ------ Generation {gen + 1}')

    # ---------- current generation
    c_rslt = rslt_data[rslt_data['Gen'] == gen]
    c_fitness = c_rslt['E_eV_atom'].to_dict()    # {ID: energy, ...}

    # ---------- load ea_data, ea_data is used only in this module
    elite_struc, elite_fitness, ea_info, ea_origin = pkl_data.load_ea_data()

    # ---------- instantiate Seclect_parents class
    logger.info('# -- select parents')
    sp = Select_parents(opt_struc_data, c_fitness, elite_struc, elite_fitness, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament()
    else:
        sp.set_roulette()

    # ---------- generate offspring by EA
    logger.info('# -- Generate structures')
    if rin.struc_mode not in ['mol', 'mol_bs']:
        _, eagen = child_gen(sp, init_struc_data, None)
    else:
        _, eagen, _ = child_gen(sp, init_struc_data, struc_mol_id)

    # ---------- select elite
    if rin.n_elite > 0:
        logger.info('# -- select elites')
        # ------ init
        all_fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}
        elite_struc = {}
        elite_fitness = {}
        # ------ Select_parents class also works for selecting elite structures
        se = Select_parents(opt_struc_data, all_fitness, None, None, rin.n_elite)
        for cid in se.ranking_dedupe:
            logger.info(f'Structure ID {cid:>6} keeps as the elite')
            elite_struc[cid] = opt_struc_data[cid]
            elite_fitness[cid] = all_fitness[cid]
    else:
        elite_struc = None
        elite_fitness = None

    # ---------- new generation
    gen += 1

    # ---------- id_queueing
    id_queueing = [i for i in range(rin.tot_struc, rin.tot_struc + rin.n_pop)]

    # ---------- ea_info
    tmp_info = pd.DataFrame(data=[[gen, rin.n_pop, rin.n_crsov, rin.n_perm,
                                   rin.n_strain, rin.n_rand, rin.n_elite,
                                   rin.crs_lat, rin.slct_func]],
                            columns=ea_info.columns)
    ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    # ------ out ea_info
    out_results.out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cid in range(rin.tot_struc, rin.tot_struc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.DataFrame([[gen, cid, eagen.operation[cid],
                                    eagen.parents[cid]]], columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ random part
    for cid in range(rin.tot_struc + rin.n_pop - rin.n_rand,
                     rin.tot_struc + rin.n_pop):
        tmp_origin = pd.DataFrame([[gen, cid, 'random', None]],
                                  columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ elite part
    if rin.n_elite > 0:
        for cid in se.ranking_dedupe:
            tmp_origin = pd.DataFrame([[gen, cid, 'elite', 'elite']],
                                      columns=ea_origin.columns)
            ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ out ea_origin
    out_results.out_ea_origin(ea_origin)

    # ---------- save ea_id_data
    ea_id_data = (gen, id_queueing, id_running)
    pkl_data.save_ea_id(ea_id_data)

    # ---------- save ea_data
    ea_data = (elite_struc, elite_fitness, ea_info, ea_origin)
    pkl_data.save_ea_data(ea_data)

    # ---------- change the value of tot_struc
    config = change_input.config_read()
    change_input.change_basic(config, 'tot_struc', rin.tot_struc + rin.n_pop)
    change_input.write_config(config)
    logger.info('# -- changed cryspy.in')
    logger.info('Changed the value of tot_struc in cryspy.in'
          f' from {rin.tot_struc} to {rin.tot_struc + rin.n_pop}')

    # ---------- status
    io_stat.set_input_common(stat, 'basic', 'tot_struc', rin.tot_struc + rin.n_pop)
    io_stat.set_common(stat, 'generation', gen)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- ext
    if rin.calc_code == 'ext':
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('out\n')
