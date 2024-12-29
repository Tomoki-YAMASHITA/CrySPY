'''
Generational change in evolutionary algorithm
'''

from logging import getLogger

import pandas as pd

from .ea_child import child_gen
from .survival import survival_fittest
from ..IO import io_stat, pkl_data
from ..IO.out_results import out_ea_info, out_ea_origin


logger = getLogger('cryspy')


def next_gen(
        rin,
        gen,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data=None,
        struc_mol_id=None,
    ):

    # ---------- log
    logger.info('# ---------- Evolutionary algorithm')
    logger.info(f'Generation {gen + 1}')

    # ---------- load ea_data. ea_data is used only in this module
    #                load hdist_data in EA-vc
    elite_struc = pkl_data.load_elite_struc()
    elite_fitness = pkl_data.load_elite_fitness()
    ea_info = pkl_data.load_ea_info()
    ea_origin = pkl_data.load_ea_origin()
    if rin.algo == 'EA-vc':
        hdist_data = pkl_data.load_hdist_data()

    # ---------- current generation
    c_rslt = rslt_data[rslt_data['Gen'] == gen]
    cgen_ids = c_rslt.index.values    # current IDs [array]

    # ---------- fitness
    if rin.algo != 'EA-vc':
        c_fitness = c_rslt['E_eV_atom'].to_dict()    # {ID: energy, ...}
        logger.debug(f'c_fitness: {c_fitness}')
    else:
        hdist = hdist_data[gen]
        c_fitness = {cid: hdist[cid] for cid in cgen_ids}
        logger.debug(f'c_fitness: {c_fitness}')
        # ------ update elite_fitness
        #        need to update elite_fitness every time hull distance is updated
        if elite_fitness is not None:
            for cid in elite_fitness:
                elite_fitness[cid] = hdist[cid]
            logger.debug(f'elite_fitness in EA-vc: {elite_fitness}')

    # ---------- survival_fittest
    logger.info('# ------ survival of the fittest')
    if rin.algo == 'EA-vc':
        # emax_ea and emin_ea are used in hdist, not in survival_fittest
        emax_ea = None
        emin_ea = None
    else:
        emax_ea = rin.emax_ea
        emin_ea = rin.emin_ea
    c_struc_data = {cid: opt_struc_data[cid] for cid in cgen_ids}
    ranking, fit_with_elite, struc_with_elite = survival_fittest(
        c_fitness,
        c_struc_data,
        elite_struc,
        elite_fitness,
        rin.n_fittest,
        rin.fit_reverse,
        emax_ea,
        emin_ea,
    )
    logger.info('ranking without duplication (including elite):')
    for cid in ranking:
            logger.info(f'Structure ID {cid:>6}, fitness: {fit_with_elite[cid]:>10.5f}')

    # ---------- generate children by EA
    logger.info('# ------ Generate children')
    pre_nstruc = len(init_struc_data)
    # init_struc_data will be updated and  saved in child_gen function
    init_struc_data, parents, operation = child_gen(
        rin,
        ranking,
        fit_with_elite,
        struc_with_elite,
        init_struc_data,
        struc_mol_id,
        nat_data,
    )

    # ---------- select elite for next generation
    if rin.n_elite > 0:
        logger.info('# ------ Select elites')
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
        # ------ ranking for all data
        ranking, _, _ = survival_fittest(
            fitness,
            opt_struc_data,
            None,
            None,
            n_elite,
            rin.fit_reverse,
            emax_ea,
            emin_ea,
        )
        if rin.algo == 'EA-vc':
            '''
            e.g.
            vert_id = [3, 8, 4] <-- assume: structure 3, 8, 4 are identical
            rin.n_elite = 2
            n_elite = len(vert_id) + rin.n_elite = 3 + 2 = 5
            ranking = [3, 6, 9, 11, 15] <-- 8 and 4 are removed because of duplication.
                                            so n_elite should be 1 + 2 = 3.
                        11 and 15 were not originally supposed to be chosen as elite

            Trim the terminals of ranking by the number of duplicated structures in vert_id
            '''
            for cid in vert_id:
                if cid not in ranking:
                    ranking.pop(-1)
        for cid in ranking:
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
    id_queueing = [i for i in range(pre_nstruc, pre_nstruc + rin.n_pop)]

    # ---------- ea_info
    if rin.algo == 'EA':
        tmp_info = pd.DataFrame(
            data=[
                [
                    gen,
                    rin.n_pop,
                    rin.n_crsov,
                    rin.n_perm,
                    rin.n_strain,
                    rin.n_rand,
                    rin.n_elite,
                    rin.crs_lat,
                    rin.slct_func
                ]
            ],
            columns=ea_info.columns
        )
    if rin.algo == 'EA-vc':
        tmp_info = pd.DataFrame(
            data=[
                [
                    gen,
                    rin.n_pop,
                    rin.n_crsov,
                    rin.n_perm,
                    rin.n_strain,
                    rin.n_add,
                    rin.n_elim,
                    rin.n_subs,
                    rin.n_rand,
                    rin.n_elite,
                    rin.crs_lat,
                    rin.slct_func,
                ]
            ],
            columns=ea_info.columns
        )
    ea_info = pd.concat([ea_info, tmp_info], axis=0, ignore_index=True)
    # ------ out ea_info
    out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    for cid in range(pre_nstruc, pre_nstruc + rin.n_pop - rin.n_rand):
        tmp_origin = pd.DataFrame([[gen, cid, operation[cid],
                                    parents[cid]]], columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ random part
    for cid in range(pre_nstruc + rin.n_pop - rin.n_rand,
                     pre_nstruc + rin.n_pop):
        tmp_origin = pd.DataFrame([[gen, cid, 'random', None]],
                                  columns=ea_origin.columns)
        ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ elite part
    if rin.n_elite > 0:
        for cid in ranking:
            tmp_origin = pd.DataFrame([[gen, cid, 'elite', 'elite']],
                                      columns=ea_origin.columns)
            ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)

    # ------ out ea_origin
    out_ea_origin(ea_origin)

    # ---------- save ea_id_data
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_gen(gen)

    # ---------- save ea_data
    pkl_data.save_elite_struc(elite_struc)
    pkl_data.save_elite_fitness(elite_fitness)
    pkl_data.save_ea_info(ea_info)
    pkl_data.save_ea_origin(ea_origin)

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_common(stat, 'generation', gen)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

