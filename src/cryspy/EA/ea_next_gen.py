'''
Generational change in evolutionary algorithm
'''

from logging import getLogger

import pandas as pd

from .calc_hull import calc_convex_hull_2d
from .ea_child import child_gen
from .gen_struc_EA.select_parents import Select_parents
from ..IO import change_input, io_stat, pkl_data
from ..IO import read_input as rin
from ..IO.out_results import out_ea_info, out_ea_origin, out_hdist

logger = getLogger('cryspy')


def next_gen(stat, init_struc_data, struc_mol_id, opt_struc_data, rslt_data, ea_id_data, ea_vc_data):
    '''
    Input:
        struc_mol_id: struc_mol_id if rin.struc_mode in ['mol', 'mol_bs'], else None
        ea_vc_data: ea_vc_data if rin.algo == 'EA-vc', else None
    '''
    # ---------- ea_id_data
    gen, id_queueing, id_running = ea_id_data

    # ---------- log
    logger.info('# ---------- Evolutionary algorithm')
    logger.info(f'# ------ Generation {gen + 1}')

    # ---------- current generation
    c_rslt = rslt_data[rslt_data['Gen'] == gen]
    if rin.algo == 'EA':
        c_fitness = c_rslt['E_eV_atom'].to_dict()    # {ID: energy, ...}

    # ---------- load ea_data. ea_data is used only in this module
    elite_struc, elite_fitness, ea_info, ea_origin = pkl_data.load_ea_data()

    # ---------- EA-vc
    if rin.algo == 'EA-vc':
        # ------ data for EA-vc
        c_ids = c_rslt.index.values    # current IDs [array]
        ef_all = rslt_data['Ef_eV_atom'].to_dict()    # formation energy of all structures
        nat_data, ratio_data, hdist_data = pkl_data.load_ea_vc_data()

        # ------ calc convex hull and hull distance
        hdist = calc_convex_hull_2d(ratio_data, ef_all, c_ids, gen)
        # -- update hdist
        out_hdist(gen, hdist, ratio_data)
        hdist_data[gen] = hdist
        ea_vc_data = (nat_data, ratio_data, hdist_data)
        pkl_data.save_ea_vc_data(ea_vc_data)

        # ------ fitness (= hull distance) of current generation
        c_fitness = {cid: hdist[cid] for cid in c_ids}
        logger.debug(f'c_fitness: {c_fitness}')

        # ------ update elite_fitness
        #        need to update elite_fitness every time hull distance (convex hull) is updated
        if elite_fitness is not None:
            for cid in elite_fitness:
                elite_fitness[cid] = hdist[cid]
            logger.debug(f'elite_fitness in EA-vc before Select_parents: {elite_fitness}')

    # ---------- instantiate Seclect_parents class
    logger.info('# -- select parents')
    sp = Select_parents(opt_struc_data, c_fitness, elite_struc, elite_fitness, rin.n_fittest)
    if rin.slct_func == 'TNM':
        sp.set_tournament()
    else:
        sp.set_roulette()

    # ---------- generate offspring by EA
    logger.info('# -- Generate structures')
    _, eagen, _ = child_gen(sp, init_struc_data, struc_mol_id, ea_vc_data)

    # ---------- select elite for next generation
    if rin.n_elite > 0:
        logger.info('# -- select elites')
        elite_struc = {}
        elite_fitness = {}
        if rin.algo == 'EA':
            fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}
            n_elite = rin.n_elite
        if rin.algo == 'EA-vc':
            '''
            In EA-vc, the num. of elite structure is
                the num. of vertices in the convex hull (hull distance < 0.001) + rin.n_elite
            '''
            vert_id = [cid for cid, value in hdist.items() if value < 0.001]
            fitness = hdist
            n_elite = len(vert_id) + rin.n_elite    # temporary
        # ------ Select_parents also works as elite selection
        se = Select_parents(opt_struc_data, fitness, None, None, n_elite)
        if rin.algo == 'EA-vc':
            '''
            e.g.
            vert_id = [3, 8, 4] <-- assume: structure 3, 8, 4 are identical
            rin.n_elite = 2
            n_elite = 3 + 2 = 5
            se.ranking_dedupe = [3, 6, 9, 11, 15] <-- 8 and 4 are removed. so n_elite should be 1 + 2 = 3.
                                11 and 15 were not originally supposed to be chosen as elite

            Trim the terminals of se.ranking_dedupe by the number of duplicated structures in vert_id
            '''
            for cid in vert_id:
                if cid not in se.ranking_dedupe:
                    se.ranking_dedupe.pop(-1)
        for cid in se.ranking_dedupe:
            logger.info(f'Structure ID {cid:>6} keeps as the elite')
            elite_struc[cid] = opt_struc_data[cid]
            elite_fitness[cid] = fitness[cid]
    else:
        elite_struc = None
        elite_fitness = None
    logger.debug(f'elite_fitness for next generation: {elite_fitness}')

    # ---------- new generation
    gen += 1

    # ---------- id_queueing
    id_queueing = [i for i in range(rin.tot_struc, rin.tot_struc + rin.n_pop)]

    # ---------- ea_info
    if rin.algo == 'EA':
        tmp_info = pd.DataFrame(data=[[gen, rin.n_pop, rin.n_crsov, rin.n_perm,
                                    rin.n_strain, rin.n_rand, rin.n_elite,
                                    rin.crs_lat, rin.slct_func]],
                                columns=ea_info.columns)
    if rin.algo == 'EA-vc':
        tmp_info = pd.DataFrame(data=[[gen, rin.n_pop,
                                    rin.n_crsov, rin.n_perm, rin.n_strain,
                                    rin.n_add, rin.n_elim, rin.n_subs,
                                    rin.n_rand, rin.n_elite,
                                    rin.crs_lat, rin.slct_func]],
                                columns=ea_info.columns)
    ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    # ------ out ea_info
    out_ea_info(ea_info)

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
    out_ea_origin(ea_origin)

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
