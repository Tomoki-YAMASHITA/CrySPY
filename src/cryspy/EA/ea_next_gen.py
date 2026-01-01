'''
Generational change in evolutionary algorithm
'''

from logging import getLogger
import os

import pandas as pd

from .ea_child import child_gen
from .survival import survival_fittest
from ..IO import io_stat, pkl_data
from ..IO.out_results import out_ea_info, out_ea_origin, out_hdist
from ..util.utility import backup_cryspy

# import later
#from ..EA.calc_hull import calc_convex_hull
#from ..util.visual_util import plot_energy_EA


logger = getLogger('cryspy')


def next_gen_EA(
        rin,
        gen,
        go_next_sg,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data=None,
        struc_mol_id=None,
        rng=None,
    ):

    # ---------- log
    logger.info(f'Done generation {gen}')

    # ---------- EA: plot
    if rin.algo == 'EA':
        from ..util.visual_util import plot_energy_EA
        # ------ plot
        fig, ax = plot_energy_EA(rslt_data, rin.ymax, rin.markersize)
        # ------ save figure
        os.makedirs('./data/energy_plot', exist_ok=True)
        gmax = rslt_data['Gen'].max()
        fname = f'./data/energy_plot/energy_EA_gen_{gmax}.{rin.fig_format}'
        fig.savefig(fname)
        logger.info(f'Energy plot for EA saved as {fname}')

    # ---------- EA-vc: calc convex hull
    if rin.algo == 'EA-vc':
        hdist_data = pkl_data.load_hdist_data()
        pd_data = pkl_data.load_pd_data()
        if gen not in hdist_data:
            logger.info(f'Calculate convex hull for generation {gen}')
            from ..EA.calc_hull import calc_convex_hull
            phase_diagram, hdist = calc_convex_hull(
                atype=rin.atype,
                gen=gen,
                end_point=rin.end_point,
                rslt_data=rslt_data,
                nat_data=nat_data,
                ymax=rin.ymax,
                show_max=rin.show_max,
                label_stable=rin.label_stable,
                vmax=rin.vmax,
                bottom_margin=rin.bottom_margin,
                markersize=rin.markersize,
                fig_format=rin.fig_format,
                emax_ea=rin.emax_ea,
                emin_ea=rin.emin_ea,
                axis_order=rin.axis_order,
            )
            out_hdist(gen, hdist, nat_data)
            hdist_data[gen] = hdist
            pd_data[gen] = phase_diagram
            pkl_data.save_hdist_data(hdist_data)
            pkl_data.save_pd_data(pd_data)

    # ---------- flag for next selection or generation
    if not go_next_sg:
        logger.info('\nEA is ready')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- check point 3
    if rin.stop_chkpt == 3:
        logger.info('\nStop at check point 3: EA is ready')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- maxgen_ea
    if 0 < rin.maxgen_ea <= gen:
        logger.info(f'\nReached maxgen_ea: {rin.maxgen_ea}')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- EA
    backup_cryspy()
    _next_gen(
        rin,
        gen,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data,
        struc_mol_id,
        rng,
    )


def _next_gen(
        rin,
        gen,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data=None,
        struc_mol_id=None,
        rng=None,
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
        c_fitness = {cid: hdist.get(cid, None) for cid in cgen_ids}
        logger.debug(f'c_fitness: {c_fitness}')
        # ------ update elite_fitness
        #        need to update elite_fitness every time hull distance is updated
        if elite_fitness is not None:
            for cid in elite_fitness:
                elite_fitness[cid] = hdist[cid]
            logger.debug(f'elite_fitness in EA-vc: {elite_fitness}')

    # ---------- natural selection
    logger.info('# ------ natural selection')
    if rin.algo == 'EA-vc':
        # emax_ea and emin_ea are used in hdist, not in survival_fittest
        emax_ea = None
        emin_ea = None
    else:
        emax_ea = rin.emax_ea
        emin_ea = rin.emin_ea
    c_struc_data = {cid: opt_struc_data[cid] for cid in cgen_ids}
    ranking, fit_with_elite, struc_with_elite = survival_fittest(
        fitness=c_fitness,
        struc_data=c_struc_data,
        elite_struc=elite_struc,
        elite_fitness=elite_fitness,
        n_fittest=rin.n_fittest,
        fit_reverse=rin.fit_reverse,
        emax_ea=emax_ea,
        emin_ea=emin_ea,
        rng=rng,
    )
    logger.info('ranking without duplication (including elite):')
    for cid in ranking:
            logger.info(f'Structure ID {cid:>6}, fitness: {fit_with_elite[cid]:>10.5f}')

    # ---------- generate children by EA
    logger.info('# ------ Generate children')
    pre_nstruc = len(init_struc_data)
    # init_struc_data will be updated and  saved in child_gen function
    init_struc_data, parents, operation = child_gen(
        rin=rin,
        ranking=ranking,
        fittest=fit_with_elite,
        struc_data=struc_with_elite,
        init_struc_data=init_struc_data,
        struc_mol_id=struc_mol_id,
        nat_data=nat_data,
        rng=rng,
    )

    # ---------- select elite for next generation
    if rin.n_elite > 0:
        logger.info('# ------ Select elites')
        elite_struc = {}
        elite_fitness = {}
        if rin.algo == 'EA':
            fitness = rslt_data['E_eV_atom'].to_dict()    # {ID: energy, ..,}
        elif rin.algo == 'EA-vc':
            fitness = hdist
        # ------ ranking for all data
        ranking, _, _ = survival_fittest(
            fitness=fitness,
            struc_data=opt_struc_data,
            elite_struc=None,
            elite_fitness=None,
            n_fittest=rin.n_elite,
            fit_reverse=rin.fit_reverse,
            emax_ea=emax_ea,
            emin_ea=emin_ea,
            rng=rng,
        )
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
    row = {
        'Gen':         gen,
        'Population':  rin.n_pop,
        'Crossover':   rin.n_crsov,
        'Permutation': rin.n_perm,
        'Strain':      rin.n_strain,
        'Random':      rin.n_rand,
        'Elite':       rin.n_elite,
        'crs_lat':     rin.crs_lat,
        'slct_func':   rin.slct_func,
    }
    if rin.algo == 'EA-vc':
        # extra fields only for EA-vc
        row.update({
            'Addition':     rin.n_add,
            'Elimination':  rin.n_elim,
            'Substitution': rin.n_subs,
        })
    ea_info.loc[len(ea_info)] = row
    # ------ out ea_info
    out_ea_info(ea_info)

    # ---------- ea_origin
    # ------ EA operation part
    rows = []
    for cid in range(pre_nstruc, pre_nstruc + rin.n_pop - rin.n_rand):
        rows.append([
            gen,
            cid,
            operation[cid],
            parents[cid],
        ])
    tmp_origin = pd.DataFrame(rows, columns=ea_origin.columns)
    ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ random part
    rows = []
    for cid in range(pre_nstruc + rin.n_pop - rin.n_rand,
                     pre_nstruc + rin.n_pop):
        rows.append([
            gen,
            cid,
            'random',
            None,
        ])
    tmp_origin = pd.DataFrame(rows, columns=ea_origin.columns)
    ea_origin = pd.concat([ea_origin, tmp_origin], axis=0, ignore_index=True)
    # ------ elite part
    if rin.n_elite > 0:
        rows = [
            [gen, cid, 'elite', 'elite']
            for cid in ranking
        ]
        tmp_origin = pd.DataFrame(rows, columns=ea_origin.columns)
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

