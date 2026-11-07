from contextlib import closing
from copy import deepcopy
from datetime import timedelta
from logging import getLogger
from logging.handlers import QueueListener
import multiprocessing as mp
import os
from queue import Empty
from time import perf_counter

import numpy as np
from pymatgen.core import Structure

from ...EA.ea_offspring import (
    EAOffspringGenerator,
    OffspringTask,
    ParentData,
    TournamentSelectionContext,
    build_roulette_selection_context,
)
from ...EA.gen_struc_EA.crossover import (
    CrossoverContext,
    CrossoverOffspringGenerator,
)
from ...EA.gen_struc_EA.permutation import (
    PermutationContext,
    PermutationOffspringGenerator,
)
from ...EA.gen_struc_EA.strain import (
    StrainContext,
    StrainOffspringGenerator,
)
from ...EA.gen_struc_EA.addition import (
    AdditionContext,
    AdditionOffspringGenerator,
    get_add_dnat_comb,
)
from ...EA.gen_struc_EA.elimination import (
    EliminationContext,
    EliminationOffspringGenerator,
    get_elim_dnat_comb,
)
from ...EA.gen_struc_EA.substitution import (
    SubstitutionContext,
    SubstitutionOffspringGenerator,
    get_subs_comb,
    apply_subs_to_nat,
)
from ...EA.natural_selection import natural_selection
from ...util.charge_neutral import prepare_cn_data
from ...util.struc_util import (
    set_mindist,
    is_comp_in_range,
    get_feasible_composition,
    precompute_feasible_N,
)
from .controller_rs import generate_structures_rs
from .worker_ea import run_worker_ea
from ..db.convert_db import raw_to_struc
from ..db.ea_db import (
    delete_ea_elites,
    insert_ea_elites,
    insert_ea_info,
    insert_ea_origin,
    select_ea_elite_ids,
    select_ea_hdist,
    select_ea_composition_bounds,
    select_ea_parent_settings,
)
from ..db.record_db import (
    Status,
    get_next_record_id,
    initialize_record,
    select_energy_results,
    select_evaluated_ids,
    select_hull_results,
    select_ids_by_generation,
    select_parent_results,
    select_status,
    select_struc,
    update_generation,
    update_init_struc,
)
from ..db.sqlite_db import connect_db


logger = getLogger('cryspy')


def plot_energy_ea(rin, generation: int):
    """Automatically plot EA energies at the specified generation interval."""

    # ---------- energy output interval
    if not (
        0 < rin.energy_output_interval
        and generation % rin.energy_output_interval == 0
    ):
        return

    import pandas as pd

    from ...util.visual_util import plot_energy_EA

    # ---------- load results
    with closing(connect_db()) as conn, conn:
        rslt_data = pd.DataFrame(
            select_energy_results(conn),
            columns=['ID', 'Gen', 'E_eV_atom'],
        ).set_index('ID')
    rslt_data['E_eV_atom'] = pd.to_numeric(rslt_data['E_eV_atom'])
    rslt_data['E_eV_atom'] = rslt_data['E_eV_atom'].replace(
        [float('inf'), -float('inf')],
        float('nan'),
    )

    # ---------- plot
    fig, ax = plot_energy_EA(rslt_data, rin.ymax, rin.markersize)

    # ---------- save figure
    os.makedirs('data/energy_plot', exist_ok=True)
    fname = f'data/energy_plot/energy_EA_gen_{generation}.{rin.fig_format}'
    fig.savefig(fname)
    logger.info(f'Energy plot for EA saved as {fname}')


def calculate_hull_ea(
    conn,
    rin,
    generation: int,
    mpl_draw: bool | None = None,
):
    """Calculate convex hull and hull distances for EA-vc."""

    from ...EA.calc_hull import calc_convex_hull

    # ---------- plot option
    if mpl_draw is None:
        mpl_draw = (
            0 < rin.energy_output_interval
            and generation % rin.energy_output_interval == 0
        )

    # ---------- output directory
    if mpl_draw:
        os.makedirs('data/convex_hull', exist_ok=True)

    # ---------- composition constraints
    comp_bounds = select_ea_composition_bounds(conn, generation)
    if comp_bounds is None:
        raise ValueError(f'Generation {generation} was not found in ea_info')
    min_comp, max_comp = comp_bounds

    # ---------- load results
    results = select_hull_results(conn, generation)
    energy_data = {cid: energy for cid, _gen, energy, _nat in results}
    nat_data = {cid: nat for cid, _gen, _energy, nat in results}

    # ---------- calculate convex hull
    logger.info(f'Calculate convex hull for generation {generation}')
    return calc_convex_hull(
        atype=rin.atype,
        gen=generation,
        ref_energies=rin.ref_energies,
        energy_data=energy_data,
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
        mpl_draw=mpl_draw,
        axis_order=rin.axis_order,
        min_comp=min_comp,
        max_comp=max_comp,
        show_comp_window=rin.show_comp_window,
    )


def load_current_generation_ids(
    conn,
    generation: int,
) -> list[int]:
    """Load structure IDs in current EA generation."""

    # ---------- current generation IDs
    current_ids = select_ids_by_generation(conn, generation)

    # ---------- check IDs
    if not current_ids:
        logger.error(
            'No structures were found for EA generation '
            f'{generation}'
        )
        raise SystemExit(1)

    # ---------- return
    return current_ids


def load_generation_results(
    conn,
    record_ids: list[int],
    hdist: dict[int, float] | None = None,
    vc: bool = False,
    struc_cache: dict[int, Structure] | None = None,
) -> tuple[
    dict[int, Structure],
    dict[int, float],
    dict[int, tuple[int, ...]] | None,
]:
    """Load optimized structures and fitness, with atom counts for EA-vc."""

    # ---------- load generation results
    results = select_parent_results(
        conn,
        record_ids,
        vc=vc,
    )
    struc_data = {}
    fitness = {}
    nat_data = {} if vc else None

    # ---------- load optimized structures
    for cid, energy, nat in results:
        # ------ reuse structures loaded during elite selection
        if struc_cache is not None and cid in struc_cache:
            struc_data[cid] = struc_cache[cid]
        else:
            # ------ load optimized structure
            raw_struc = select_struc(conn, cid, initial=False)
            if raw_struc is None:
                continue
            if raw_struc['lattice'] is None:
                continue
            if raw_struc['frac_coords'] is None:
                continue

            struc_data[cid] = raw_to_struc(raw_struc)
        fitness[cid] = energy if hdist is None else hdist.get(cid, np.nan)
        if nat_data is not None:
            nat_data[cid] = nat

    # ---------- check generation results
    missing_ids = [
        cid
        for cid in sorted(set(record_ids) - set(struc_data))
        if select_status(conn, cid) not in (
            Status.ERROR,
            Status.MINDIST,
            Status.OMIT,
        )
    ]
    if missing_ids:
        logger.error(
            'Failed to load generation results: '
            f'missing IDs = {missing_ids}'
        )
        raise SystemExit(1)

    # ---------- return
    return struc_data, fitness, nat_data


def filter_parent_pool_by_nat(
    rin,
    struc_data: dict[int, Structure],
    fitness: dict[int, float],
    nat_data: dict[int, tuple[int, ...]],
):
    """Filter parent candidates by atom-count and composition bounds."""

    # ---------- initialize
    filtered_struc = {}
    filtered_fitness = {}

    # ---------- filter candidates
    for cid, struc in struc_data.items():
        nat = nat_data.get(cid)
        if nat is None or sum(nat) <= 0:
            continue
        if not all(
            ll <= n_i <= ul
            for n_i, ll, ul in zip(nat, rin.ll_nat, rin.ul_nat)
        ):
            continue
        # ------ composition constraints
        if not is_comp_in_range(nat, rin.min_comp, rin.max_comp):
            continue

        # ------ valid fitness
        value = fitness.get(cid)
        if value is None or not np.isfinite(value):
            continue

        filtered_struc[cid] = struc
        filtered_fitness[cid] = value

    # ---------- return
    return filtered_struc, filtered_fitness


def prepare_parent_pool(
    conn,
    rin,
    generation: int,
    rng: np.random.Generator,
):
    """Prepare parent pool for EA and EA-vc."""

    # ---------- hull distance
    hdist = select_ea_hdist(conn, generation) if rin.algo == 'EA-vc' else None

    # ---------- current generation
    current_ids = load_current_generation_ids(conn, generation)
    struc_data, fitness, nat_data = load_generation_results(
        conn,
        current_ids,
        hdist=hdist,
        vc=rin.algo == 'EA-vc',
    )
    if rin.algo == 'EA-vc':
        filtered_struc_data, filtered_fitness = filter_parent_pool_by_nat(
            rin,
            struc_data,
            fitness,
            nat_data,
        )
    else:
        filtered_struc_data = struc_data
        filtered_fitness = fitness

    # ---------- elite
    elite_ids = select_ea_elite_ids(conn, generation)
    reselect_reasons = []
    if any(
        select_status(conn, cid) == Status.OMIT
        for cid in elite_ids
    ):
        reselect_reasons.append('saved elite was omitted')

    if rin.algo == 'EA-vc':
        if rin.n_elite == 0:
            elite_ids = []
            reselect_reasons.clear()
            save_ea_elite_ids(conn, generation, elite_ids)
        elif generation > 1:
            previous_settings = select_ea_parent_settings(conn, generation)
            current_settings = (
                rin.ll_nat,
                rin.ul_nat,
                rin.n_elite,
                rin.min_comp,
                rin.max_comp,
            )
            if previous_settings != current_settings:
                reselect_reasons.append(
                    'atom-count bounds, composition bounds or n_elite changed'
                )
            if not elite_ids:
                reselect_reasons.append('no saved elites')
            elif any(
                not np.isfinite(hdist.get(cid, np.nan))
                for cid in elite_ids
            ):
                reselect_reasons.append('saved elite has no valid hull distance')

    # Structures are cached only when elites are reselected.
    elite_struc_cache = None
    if reselect_reasons:
        # ------ reselect elites from past generations
        logger.info(
            f'Reselect elites from generations <= {generation - 1}: '
            + ', '.join(reselect_reasons)
        )
        past_ids = select_evaluated_ids(conn, generation - 1)
        # Use past candidates with hull distances from the current generation.
        past_struc, past_fitness = prepare_elite_pool(
            conn,
            rin,
            generation,
            record_ids=past_ids,
        )

        # ------ structure loader
        struc_loader = prepare_struc_loader(conn)

        # ------ elite RNG
        # Reselection must not advance the parent-selection RNG on restart.
        elite_rng = deepcopy(rng)

        # ------ select elite
        elite_ids, elite_struc_cache = select_elite_ids(
            rin,
            past_struc,
            past_fitness,
            elite_rng,
            struc_loader=struc_loader,
        )

        # ------ update latest elite record
        save_ea_elite_ids(conn, generation, elite_ids)

    if elite_ids:
        elite_struc, elite_fitness, elite_nat_data = load_generation_results(
            conn,
            elite_ids,
            hdist=hdist,
            vc=rin.algo == 'EA-vc',
            struc_cache=elite_struc_cache,
        )
        if rin.algo == 'EA-vc':
            filtered_elite_struc, filtered_elite_fitness = filter_parent_pool_by_nat(
                rin,
                elite_struc,
                elite_fitness,
                elite_nat_data,
            )
            nat_data.update({
                cid: elite_nat_data[cid]
                for cid in filtered_elite_struc
            })
        else:
            filtered_elite_struc = elite_struc
            filtered_elite_fitness = elite_fitness
    else:
        filtered_elite_struc = None
        filtered_elite_fitness = None

    # ---------- log valid parent candidates
    if rin.algo == 'EA-vc':
        n_elite_candidates = (
            0 if filtered_elite_struc is None else len(filtered_elite_struc)
        )
        logger.info(
            f'Current parent candidates: {len(filtered_struc_data)}'
        )
        logger.info(
            f'Elite parent candidates: {n_elite_candidates}'
        )

    # ---------- natural selection
    logger.info('# ------ natural selection')
    ranking, fit_with_elite, struc_with_elite = natural_selection(
        fitness=filtered_fitness,
        struc_data=filtered_struc_data,
        elite_struc=filtered_elite_struc,
        elite_fitness=filtered_elite_fitness,
        n_fittest=rin.n_fittest,
        fit_reverse=rin.fit_reverse,
        emax_ea=None if rin.algo == 'EA-vc' else rin.emax_ea,
        emin_ea=None if rin.algo == 'EA-vc' else rin.emin_ea,
        rng=rng,
    )

    # ---------- log parent ranking
    logger.info('ranking without duplication (including elite):')
    for cid in ranking:
        logger.info(
            f'Structure ID {cid:>6}, fitness: {fit_with_elite[cid]:>10.5f}'
        )

    # ---------- check parent pool
    min_required = 2 if rin.n_crsov > 0 else 1
    if len(ranking) < min_required:
        logger.error(
            'Not enough parent candidates after natural selection: '
            f'required {min_required}, got {len(ranking)}'
        )
        if rin.algo == 'EA-vc':
            logger.error(
                'Parent candidates after nat/composition filtering: '
                f'current={len(filtered_struc_data)}, elite={n_elite_candidates}; '
                f'll_nat={rin.ll_nat}, ul_nat={rin.ul_nat}; '
                f'min_comp={rin.min_comp}, max_comp={rin.max_comp}'
            )
        raise SystemExit(1)

    # ---------- return
    selected_struc_data = {cid: struc_with_elite[cid] for cid in ranking}
    selected_nat_data = None
    if nat_data is not None:
        selected_nat_data = {cid: nat_data[cid] for cid in ranking}
    return ranking, fit_with_elite, selected_struc_data, selected_nat_data


def prepare_parent_data(
    struc_data: dict[int, Structure],
) -> ParentData:
    """Prepare parent structure data."""

    # ---------- parent data
    return ParentData(
        structures=struc_data,
    )


def prepare_selection_context(
    rin,
    ranking: list[int],
    fitness: dict[int, float],
):
    """Prepare parent selection context."""

    # ---------- tournament
    if rin.slct_func == 'TNM':
        return TournamentSelectionContext(
            parent_ids=tuple(ranking),
            tournament_size=rin.t_size,
        )

    return build_roulette_selection_context(
        parent_ids=tuple(ranking),
        fitness=fitness,
        a_rlt=rin.a_rlt,
        b_rlt=rin.b_rlt,
        fit_reverse=rin.fit_reverse,
    )


def prepare_crossover_context(
    rin,
    mindist,
    feasible_N: list | None = None,
    cn_data: dict | None = None,
) -> CrossoverContext:
    """Prepare crossover context."""

    # ---------- crossover context
    return CrossoverContext(
        atype=rin.atype,
        nat=rin.nat,
        mindist=mindist,
        crs_lat=rin.crs_lat,
        nat_diff_tole=rin.nat_diff_tole,
        maxcnt_ea=rin.maxcnt_ea,
        vc=rin.algo == 'EA-vc',
        ll_nat=rin.ll_nat,
        ul_nat=rin.ul_nat,
        cn_comb=cn_data['cn_comb'] if cn_data is not None else None,
        feasible_N=feasible_N,
        charge=rin.charge if rin.algo == 'EA-vc' else None,
        cn_data=cn_data,
        min_comp=rin.min_comp if rin.algo == 'EA-vc' else None,
        max_comp=rin.max_comp if rin.algo == 'EA-vc' else None,
    )


def prepare_permutation_parent_ids(
    ranking: list[int],
    nat_data: dict[int, tuple[int, ...]],
) -> list[int]:
    """Prepare EA-vc parent candidates for permutation."""

    # ---------- permutation candidates
    parent_ids = []
    for pid in ranking:
        nat = nat_data[pid]
        if sum(n > 0 for n in nat) >= 2:
            parent_ids.append(pid)

    # ---------- check candidates
    if not parent_ids:
        logger.error('No feasible parent candidates for permutation')
        raise SystemExit(1)

    # ---------- return
    return parent_ids


def prepare_permutation_context(
    rin,
    mindist,
) -> PermutationContext:
    """Prepare permutation context."""

    # ---------- permutation context
    return PermutationContext(
        atype=rin.atype,
        mindist=mindist,
        ntimes=rin.ntimes,
        maxcnt_ea=rin.maxcnt_ea,
    )


def prepare_strain_context(
    rin,
    mindist,
) -> StrainContext:
    """Prepare strain context."""

    # ---------- strain context
    return StrainContext(
        atype=rin.atype,
        mindist=mindist,
        sigma_st=rin.sigma_st,
        maxcnt_ea=rin.maxcnt_ea,
    )


def prepare_addition_context(
    rin,
    mindist,
    ranking: list[int],
    nat_data: dict[int, tuple[int, ...]],
    cn_data: dict | None = None,
) -> AdditionContext:
    """Prepare addition context under atom-count and composition bounds."""

    # ---------- addition candidates
    add_dnat_map = {}
    for pid in ranking:
        parent_nat = nat_data[pid]
        dnat_comb = get_add_dnat_comb(
            rin.ul_nat,
            rin.add_max,
            parent_nat,
            cn_comb=cn_data['cn_comb'] if cn_data is not None else None,
            charge=rin.charge,
        )
        # ------ composition constraints
        dnat_ok = []
        for dnat in dnat_comb:
            new_nat = tuple(a + b for a, b in zip(parent_nat, dnat))
            if is_comp_in_range(new_nat, rin.min_comp, rin.max_comp):
                dnat_ok.append(dnat)
        if dnat_ok:
            add_dnat_map[pid] = dnat_ok

    # ---------- check candidates
    if not add_dnat_map:
        logger.error(
            'No feasible parent candidates for addition under current '
            'nat/composition constraints: '
            f'ul_nat={rin.ul_nat}, add_max={rin.add_max}, '
            f'min_comp={rin.min_comp}, max_comp={rin.max_comp}'
        )
        raise SystemExit(1)

    # ---------- addition context
    return AdditionContext(
        atype=rin.atype,
        mindist=mindist,
        add_dnat_map=add_dnat_map,
        maxcnt_ea=rin.maxcnt_ea,
    )


def prepare_elimination_context(
    rin,
    ranking: list[int],
    nat_data: dict[int, tuple[int, ...]],
    cn_data: dict | None = None,
) -> EliminationContext:
    """Prepare elimination context under atom-count and composition bounds."""

    # ---------- elimination candidates
    elim_dnat_map = {}
    for pid in ranking:
        parent_nat = nat_data[pid]
        dnat_comb = get_elim_dnat_comb(
            rin.ll_nat,
            rin.elim_max,
            parent_nat,
            cn_comb=cn_data['cn_comb'] if cn_data is not None else None,
            charge=rin.charge,
        )
        # ------ composition constraints
        dnat_ok = []
        for dnat in dnat_comb:
            new_nat = tuple(a - b for a, b in zip(parent_nat, dnat))
            if sum(new_nat) > 0 and is_comp_in_range(
                new_nat, rin.min_comp, rin.max_comp
            ):
                dnat_ok.append(dnat)
        if dnat_ok:
            elim_dnat_map[pid] = dnat_ok

    # ---------- check candidates
    if not elim_dnat_map:
        logger.error(
            'No feasible parent candidates for elimination under current '
            'nat/composition constraints: '
            f'll_nat={rin.ll_nat}, elim_max={rin.elim_max}, '
            f'min_comp={rin.min_comp}, max_comp={rin.max_comp}'
        )
        raise SystemExit(1)

    # ---------- elimination context
    return EliminationContext(
        atype=rin.atype,
        elim_dnat_map=elim_dnat_map,
    )


def prepare_substitution_context(
    rin,
    mindist,
    ranking: list[int],
    nat_data: dict[int, tuple[int, ...]],
    cn_data: dict | None = None,
) -> SubstitutionContext:
    """Prepare substitution context under atom-count and composition bounds."""

    # ---------- substitution candidates
    subs_comb_map = {}
    cn_comb = cn_data['cn_comb'] if cn_data is not None else None
    cn_set = {tuple(row) for row in cn_comb} if cn_comb is not None else None
    for pid in ranking:
        parent_nat = nat_data[pid]
        subs_comb = get_subs_comb(
            rin.atype,
            rin.ll_nat,
            rin.ul_nat,
            rin.subs_max,
            parent_nat,
            cn_set=cn_set,
            charge=rin.charge,
        )
        # ------ composition constraints
        subs_ok = []
        for subs_list in subs_comb:
            new_nat = apply_subs_to_nat(parent_nat, subs_list, rin.atype)
            if new_nat is not None and is_comp_in_range(
                new_nat, rin.min_comp, rin.max_comp
            ):
                subs_ok.append(subs_list)
        if subs_ok:
            subs_comb_map[pid] = subs_ok

    # ---------- check candidates
    if not subs_comb_map:
        logger.error(
            'No feasible parent candidates for substitution under current '
            'nat/composition constraints: '
            f'll_nat={rin.ll_nat}, ul_nat={rin.ul_nat}, '
            f'subs_max={rin.subs_max}, '
            f'min_comp={rin.min_comp}, max_comp={rin.max_comp}'
        )
        raise SystemExit(1)

    # ---------- substitution context
    return SubstitutionContext(
        atype=rin.atype,
        mindist=mindist,
        subs_comb_map=subs_comb_map,
        maxcnt_ea=rin.maxcnt_ea,
    )


def create_offspring_tasks(
    conn,
    rin,
    generation: int,
) -> list[OffspringTask]:
    """Create offspring generation tasks."""

    # ---------- initialize
    tasks = []
    cid = get_next_record_id(conn)

    # ---------- operation counts
    operation_counts = (
        ('crossover', rin.n_crsov),
        ('permutation', rin.n_perm),
        ('strain', rin.n_strain),
    )
    if rin.algo == 'EA-vc':
        operation_counts += (
            ('addition', rin.n_add),
            ('elimination', rin.n_elim),
            ('substitution', rin.n_subs),
        )

    # ---------- create tasks
    for operation, count in operation_counts:
        for _ in range(count):
            initialize_record(
                conn,
                cid,
                Status.GENERATING,
                generation=generation,
            )
            tasks.append(
                OffspringTask(
                    cid=cid,
                    generation=generation,
                    operation=operation,
                )
            )
            cid += 1

    # ---------- return
    return tasks


def prepare_random_records_ea(
    conn,
    random_ids: list[int],
    generation: int,
) -> list[int]:
    """Prepare records for EA random structure generation."""

    # ---------- initialize
    generate_ids = []

    # ---------- prepare records
    for cid in random_ids:
        status = select_status(conn, cid)

        # ------ existing record
        if status is not None:
            if status == Status.GENERATING:
                generate_ids.append(cid)
            continue

        # ------ initialize missing record
        initialize_record(
            conn,
            cid,
            Status.GENERATING,
            generation=generation,
        )
        generate_ids.append(cid)

    # ---------- return
    return generate_ids


def prepare_offspring_generator(
    rin,
    parent_data: ParentData,
    selection_context,
    selection_context_pm,
    crossover_context: CrossoverContext,
    permutation_context: PermutationContext,
    strain_context: StrainContext,
    selection_context_ad,
    addition_context: AdditionContext | None,
    selection_context_el,
    elimination_context: EliminationContext | None,
    selection_context_sb,
    substitution_context: SubstitutionContext | None,
) -> EAOffspringGenerator:
    """Prepare EA offspring generator."""

    # ---------- operation generators
    operation_generators = {
        'crossover': CrossoverOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context,
            context=crossover_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate,
        'permutation': PermutationOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context_pm,
            context=permutation_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate,
        'strain': StrainOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context,
            context=strain_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate,
    }

    # ---------- addition
    if addition_context is not None:
        operation_generators['addition'] = AdditionOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context_ad,
            context=addition_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate

    # ---------- elimination
    if elimination_context is not None:
        operation_generators['elimination'] = EliminationOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context_el,
            context=elimination_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate

    # ---------- substitution
    if substitution_context is not None:
        operation_generators['substitution'] = SubstitutionOffspringGenerator(
            parent_data=parent_data,
            selection_context=selection_context_sb,
            context=substitution_context,
            max_parent_attempts=rin.max_parent_attempts,
        ).generate

    # ---------- offspring generator
    return EAOffspringGenerator(
        operation_generators=operation_generators,
    )


def run_controller_ea(
    tasks: list[OffspringTask],
    task_queue,
    result_queue,
    stop_event,
    processes: list[mp.Process],
):
    """Control EA offspring generation."""

    # ---------- initialize
    task_iter = iter(tasks)
    num_active = 0
    worker_failed = False
    stop_requested = False
    results = []

    # ---------- initial tasks
    for _ in processes:
        if os.path.isfile('STOP_CRYSPY_HT'):
            stop_requested = True
            logger.info('Stop file detected: STOP_CRYSPY_HT')
            break

        try:
            task = next(task_iter)
        except StopIteration:
            break
        task_queue.put(task)
        num_active += 1

    # ---------- process results
    while num_active > 0 or worker_failed:
        result = None
        try:
            # ------ get result
            result = result_queue.get(timeout=1.0)

        except Empty:
            pass

        # ------ check workers
        if not worker_failed:
            stopped_processes = [
                process for process in processes
                if process.exitcode is not None
            ]
            if stopped_processes:
                for process in stopped_processes:
                    logger.error(
                        f'{process.name} stopped: '
                        f'exitcode = {process.exitcode}'
                    )

                # ------ drain workers
                worker_failed = True
                stop_event.set()
                for _ in processes:
                    task_queue.put(None)

        # ------ no result
        if result is None:
            if (
                worker_failed
                and all(
                    process.exitcode is not None
                    for process in processes
                )
            ):
                break
            continue

        task, offspring_result, error_message = result
        num_active -= 1

        # ------ check generation error
        if error_message is not None:
            raise RuntimeError(
                f'Offspring generation failed: '
                f'ID {task.cid}: {error_message}'
            )

        # ------ store result
        results.append((task, offspring_result))

        # ------ stop assigning tasks
        if worker_failed:
            continue

        # ------ check stop file
        if not stop_requested and os.path.isfile('STOP_CRYSPY_HT'):
            stop_requested = True
            logger.info('Stop file detected: STOP_CRYSPY_HT')
        if stop_requested:
            continue

        # ------ queue next task
        try:
            task = next(task_iter)
        except StopIteration:
            continue
        task_queue.put(task)
        num_active += 1

    # ---------- worker failure
    if worker_failed:
        raise RuntimeError(
            'Worker stopped before returning all offspring'
        ) from None

    # ---------- return
    return results, stop_requested


def save_offspring_results(
    conn,
    rin,
    results,
) -> None:
    """Save EA offspring generation results."""

    # ---------- save results
    for task, offspring_result in results:
        parent_ids = offspring_result.parent_ids
        parent1_id = parent_ids[0] if len(parent_ids) > 0 else None
        parent2_id = parent_ids[1] if len(parent_ids) > 1 else None

        update_init_struc(
            conn,
            task.cid,
            offspring_result.structure,
            rin.atype,
            rin.symprec,
            Status.WAITING,
            vc=rin.algo == 'EA-vc',
        )
        insert_ea_origin(
            conn,
            task.cid,
            task.generation,
            task.operation,
            parent1_id,
            parent2_id,
            offspring_result.parent_attempt_count,
        )


def save_ea_generation_info(
    conn,
    rin,
    generation: int,
) -> None:
    """Save EA generation information."""

    # ---------- save generation information
    insert_ea_info(
        conn,
        generation,
        rin.n_pop,
        rin.n_crsov,
        rin.n_perm,
        rin.n_strain,
        rin.n_rand,
        rin.n_elite,
        rin.crs_lat,
        rin.slct_func,
        n_addition=rin.n_add if rin.algo == 'EA-vc' else None,
        n_elimination=rin.n_elim if rin.algo == 'EA-vc' else None,
        n_substitution=rin.n_subs if rin.algo == 'EA-vc' else None,
        min_comp=rin.min_comp if rin.algo == 'EA-vc' else None,
        max_comp=rin.max_comp if rin.algo == 'EA-vc' else None,
        ll_nat=rin.ll_nat if rin.algo == 'EA-vc' else None,
        ul_nat=rin.ul_nat if rin.algo == 'EA-vc' else None,
    )


def save_ea_elite_ids(
    conn,
    generation: int,
    elite_ids: list[int],
) -> None:
    """Save EA elite IDs."""

    # ---------- replace elite IDs
    # An empty list clears the elite record for this generation.
    delete_ea_elites(conn, generation)
    insert_ea_elites(
        conn,
        generation,
        elite_ids,
    )


def prepare_elite_pool(
    conn,
    rin,
    generation: int,
    record_ids: list[int] | None = None,
):
    """Prepare elite candidates without loading structures."""

    # ---------- hull distances
    hdist = select_ea_hdist(conn, generation) if rin.algo == 'EA-vc' else None

    # ---------- load candidates
    # Use all evaluated candidates unless candidate IDs are explicitly supplied.
    if record_ids is None:
        if rin.algo == 'EA-vc':
            record_ids = sorted(hdist)
        else:
            record_ids = select_evaluated_ids(conn, generation)

    results = select_parent_results(
        conn,
        record_ids,
        vc=rin.algo == 'EA-vc',
    )

    # ---------- check candidate results
    result_ids = {cid for cid, _energy, _nat in results}
    missing_ids = [
        cid
        for cid in sorted(set(record_ids) - result_ids)
        if select_status(conn, cid) not in (
            Status.ERROR,
            Status.MINDIST,
            Status.OMIT,
        )
    ]
    if missing_ids:
        logger.error(
            'Failed to load generation results: '
            f'missing IDs = {missing_ids}'
        )
        raise SystemExit(1)

    # ---------- filter candidates
    fitness = {}
    for cid, energy, nat in results:
        if rin.algo == 'EA-vc':
            if nat is None or sum(nat) <= 0:
                continue
            if not all(
                ll <= n_i <= ul
                for n_i, ll, ul in zip(nat, rin.ll_nat, rin.ul_nat)
            ):
                continue
            # ------ composition constraints
            if not is_comp_in_range(nat, rin.min_comp, rin.max_comp):
                continue

        value = energy if hdist is None else hdist.get(cid, np.nan)
        if value is None or not np.isfinite(value):
            continue
        fitness[cid] = value

    # ---------- log valid elite candidates
    if rin.algo == 'EA-vc':
        logger.info(
            f'Elite selection candidates: {len(fitness)}'
        )

    # ---------- return
    # Keep the structure pool empty to avoid loading all candidate structures.
    # natural_selection loads only the structures needed via struc_loader.
    return {}, fitness


def prepare_struc_loader(conn):
    """Prepare a loader for optimized structures."""

    # ---------- structure loader
    def struc_loader(cid: int) -> Structure:
        # ------ load optimized structure
        raw_struc = select_struc(conn, cid, initial=False)

        # ------ check optimized structure
        if (
            raw_struc is None
            or raw_struc['lattice'] is None
            or raw_struc['frac_coords'] is None
        ):
            logger.error(
                f'Failed to load optimized structure: ID = {cid}'
            )
            raise SystemExit(1)

        # ------ convert structure
        return raw_to_struc(raw_struc)

    # ---------- return
    # Return the function bound to conn; no structures are loaded yet.
    return struc_loader


def select_elite_ids(
    rin,
    struc_data: dict[int, Structure],
    fitness: dict[int, float],
    rng: np.random.Generator,
    struc_loader=None,
) -> tuple[list[int], dict[int, Structure]]:
    """Select elite IDs and structures."""

    # ---------- no elite
    if rin.n_elite <= 0:
        return [], {}

    # ---------- select elite
    logger.info('# ------ Select elites')
    ranking, _, struc_with_elite = natural_selection(
        fitness=fitness,
        struc_data=struc_data,
        elite_struc=None,
        elite_fitness=None,
        n_fittest=rin.n_elite,
        fit_reverse=rin.fit_reverse,
        emax_ea=None if rin.algo == 'EA-vc' else rin.emax_ea,
        emin_ea=None if rin.algo == 'EA-vc' else rin.emin_ea,
        rng=rng,
        struc_loader=struc_loader,
    )

    # ---------- log elite
    for cid in ranking:
        logger.info(f'Structure ID {cid:>6} keeps as the elite')

    # ---------- return
    # Keep selected structures for reuse when loading elite results.
    elite_struc = {cid: struc_with_elite[cid] for cid in ranking}
    return ranking, elite_struc


def generate_offspring_ea(
    rin,
    tasks: list[OffspringTask],
    offspring_generator: EAOffspringGenerator,
):
    """Generate EA offspring with workers."""

    # ---------- nothing to generate
    if not tasks:
        return [], False

    # ---------- check stop file
    if os.path.isfile('STOP_CRYSPY_HT'):
        logger.info('Stop file detected: STOP_CRYSPY_HT')
        return [], True

    # ---------- log
    logger.info('# ------ Generate children')

    # ---------- number of workers
    njob_struc = rin.njob if rin.njob_struc is None else rin.njob_struc
    num_workers = min(njob_struc, len(tasks))
    logger.info(f'Number of EA generation workers: {num_workers}')

    # ---------- queues and event
    task_queue = mp.Queue()
    result_queue = mp.Queue()
    log_queue = mp.Queue()
    stop_event = mp.Event()

    # ---------- logging
    listener = QueueListener(
        log_queue,
        *logger.handlers,
        respect_handler_level=True,
    )
    listener_started = False

    # ---------- initialize workers
    processes = []

    try:
        # ---------- start workers
        for worker_id in range(num_workers):
            process = mp.Process(
                target=run_worker_ea,
                args=(
                    offspring_generator,
                    rin.seed,
                    task_queue,
                    result_queue,
                    stop_event,
                    log_queue,
                    logger.level,
                ),
                name=f'worker-ea-{worker_id + 1}',
            )
            process.start()
            processes.append(process)
            logger.info(
                f'Launched {process.name}: PID = {process.pid}'
            )

        # ---------- start listener
        listener.start()
        listener_started = True

        # ---------- run controller
        results, stop_requested = run_controller_ea(
            tasks,
            task_queue,
            result_queue,
            stop_event,
            processes,
        )

        # ---------- stop workers
        for _ in processes:
            task_queue.put(None)

        # ---------- wait workers
        for process in processes:
            process.join()

    except BaseException as error:
        # ---------- cancel task queue join
        task_queue.cancel_join_thread()

        # ---------- terminate workers
        for process in processes:
            if process.is_alive():
                process.terminate()
        for process in processes:
            process.join()

        # ---------- execution error
        if isinstance(error, Exception):
            logger.exception(
                'EA offspring generation failed'
            )
            raise SystemExit(1)

        raise

    finally:
        # ---------- stop logging
        if listener_started:
            listener.stop()

        # ---------- close queues
        task_queue.close()
        result_queue.close()
        log_queue.close()
        task_queue.join_thread()
        result_queue.join_thread()
        log_queue.join_thread()

    # ---------- check exit codes
    failed_processes = [
        process for process in processes
        if process.exitcode != 0
    ]
    if failed_processes:
        for process in failed_processes:
            logger.error(
                f'{process.name} failed: '
                f'exitcode = {process.exitcode}'
            )
        raise SystemExit(1)

    # ---------- return
    return results, stop_requested


def generate_and_save_offspring_ea(
    conn,
    rin,
    tasks: list[OffspringTask],
    offspring_generator: EAOffspringGenerator,
    elite_pool_struc: dict[int, Structure],
    elite_pool_fitness: dict[int, float],
    rng: np.random.Generator,
) -> list[int]:
    """Generate and save EA offspring."""

    # ---------- generate offspring
    results, stop_requested = generate_offspring_ea(
        rin,
        tasks,
        offspring_generator,
    )

    # ---------- select elite
    elite_ids = []
    if not stop_requested:
        # ------ structure loader
        struc_loader = prepare_struc_loader(conn)
        elite_ids, _ = select_elite_ids(
            rin,
            elite_pool_struc,
            elite_pool_fitness,
            rng,
            struc_loader=struc_loader,
        )

    try:
        # ---------- save results
        save_offspring_results(
            conn,
            rin,
            results,
        )

        conn.commit()

    except Exception:
        # ---------- rollback
        conn.rollback()
        raise

    # ---------- stop
    if stop_requested:
        # Leave this generation incomplete for the existing restart recovery.
        raise SystemExit()

    # ---------- return
    return elite_ids


def initialize_ea(
    rin,
    cids: list[int],
) -> None:
    """Initialize fixed-composition EA."""

    # ---------- generation
    generation = 1
    logger.info('# ---------- Initialize evolutionary algorithm')
    logger.info(f'# ------ Generation {generation}')

    # ---------- time
    time_start = perf_counter()

    # ---------- initialize database connection
    conn = None
    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- prepare random records
        generate_ids = prepare_random_records_ea(
            conn,
            cids,
            generation,
        )
        conn.commit()

        # ---------- generate random structures
        generate_structures_rs(
            rin,
            generate_ids,
        )


        # ---------- update generation
        update_generation(
            conn,
            cids,
            generation,
        )

        # ---------- check generation IDs
        generation_ids = set(
            select_ids_by_generation(
                conn,
                generation,
            )
        )
        expected_ids = set(cids)
        if generation_ids != expected_ids:
            missing_ids = sorted(expected_ids - generation_ids)
            unexpected_ids = sorted(generation_ids - expected_ids)
            logger.error(
                f'Initial EA population does not match n_pop: '
                f'missing IDs = {missing_ids}, '
                f'unexpected IDs = {unexpected_ids}'
            )
            raise SystemExit(1)

        # ---------- save random origins
        for cid in cids:
            insert_ea_origin(
                conn,
                cid,
                generation,
                'random',
            )

        # ---------- save generation information
        insert_ea_info(
            conn,
            generation,
            rin.n_pop,
            0,
            0,
            0,
            rin.n_pop,
            0,
            rin.crs_lat,
            rin.slct_func,
            n_addition=0 if rin.algo == 'EA-vc' else None,
            n_elimination=0 if rin.algo == 'EA-vc' else None,
            n_substitution=0 if rin.algo == 'EA-vc' else None,
            min_comp=rin.min_comp if rin.algo == 'EA-vc' else None,
            max_comp=rin.max_comp if rin.algo == 'EA-vc' else None,
            ll_nat=rin.ll_nat if rin.algo == 'EA-vc' else None,
            ul_nat=rin.ul_nat if rin.algo == 'EA-vc' else None,
        )

        conn.commit()
    except Exception:
        if conn is not None:
            conn.rollback()
        raise
    finally:
        if conn is not None:
            conn.close()

    # ---------- time
    etime = timedelta(seconds=perf_counter() - time_start)
    logger.info(
        f'Elapsed time for generation {generation} structure generation: {etime}'
    )


def generate_next_generation_ea(
    rin,
    current_generation: int,
) -> None:
    """Generate and save the next EA generation."""

    # ---------- generation
    next_generation = current_generation + 1
    logger.info('# ---------- Evolutionary algorithm')
    logger.info(f'# ------ Generation {next_generation}')

    # ---------- time
    time_start = perf_counter()

    # ---------- initialize database connection
    conn = None

    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- controller RNG
        if rin.seed is None:
            rng = np.random.default_rng()
        else:
            seed_sequence = np.random.SeedSequence(
                [
                    rin.seed,
                    current_generation,
                ]
            )
            rng = np.random.default_rng(seed_sequence)

        # ---------- charge neutrality
        cn_data = None
        if rin.algo == 'EA-vc' and rin.charge is not None:
            cn_data = prepare_cn_data(
                rin.ll_nat,
                rin.ul_nat,
                rin.charge,
                cn_mode=rin.cn_mode,
                max_cn_grid_points=rin.max_cn_grid_points,
                min_comp=rin.min_comp,
                max_comp=rin.max_comp,
            )

        # ---------- parent pool
        # selected_nat_data is used only for EA-vc; None for EA.
        ranking, fitness, selected_struc_data, selected_nat_data = prepare_parent_pool(
            conn,
            rin,
            current_generation,
            rng,
        )

        # ---------- mindist
        logger.info('# ------ mindist')
        mindist = set_mindist(
            rin.atype,
            rin.mindist,
            rin.mindist_factor,
            rin.struc_mode,
        )

        # ---------- log eligible candidates
        logger.info(
            f'eligible candidates for crossover and strain: {len(ranking)}'
        )

        # ---------- composition constraints
        feasible_N = None
        if rin.algo == 'EA-vc' and (
            rin.min_comp is not None or rin.max_comp is not None
        ):
            feasible_comp = get_feasible_composition(rin.min_comp, rin.max_comp)
            feasible_N = precompute_feasible_N(rin.ll_nat, rin.ul_nat, feasible_comp)
            logger.info(
                'Feasible total atom counts for composition constraints: '
                f'{len(feasible_N)}'
            )

        # ---------- generation contexts
        parent_data = prepare_parent_data(selected_struc_data)
        selection_context = prepare_selection_context(
            rin,
            ranking,
            fitness,
        )
        crossover_context = prepare_crossover_context(
            rin,
            mindist,
            feasible_N=feasible_N,
            cn_data=cn_data,
        )
        permutation_context = prepare_permutation_context(
            rin,
            mindist,
        )
        strain_context = prepare_strain_context(
            rin,
            mindist,
        )

        # ---------- selection context for permutation
        selection_context_pm = selection_context
        if rin.algo == 'EA-vc' and rin.n_perm > 0:
            perm_parent_ids = prepare_permutation_parent_ids(
                ranking,
                selected_nat_data,
            )
            selection_context_pm = prepare_selection_context(
                rin,
                perm_parent_ids,
                fitness,
            )
        if rin.n_perm > 0:
            logger.info(
                'eligible candidates for permutation: '
                f'{len(selection_context_pm.parent_ids)}'
            )

        # ---------- addition contexts for EA-vc
        addition_context = None
        selection_context_ad = None
        if rin.algo == 'EA-vc' and rin.n_add > 0:
            addition_context = prepare_addition_context(
                rin,
                mindist,
                ranking,
                selected_nat_data,
                cn_data=cn_data,
            )
            add_parent_ids = list(addition_context.add_dnat_map)
            logger.info(
                f'eligible candidates for addition: {len(add_parent_ids)}'
            )
            selection_context_ad = prepare_selection_context(
                rin,
                add_parent_ids,
                fitness,
            )

        # ---------- elimination contexts for EA-vc
        elimination_context = None
        selection_context_el = None
        if rin.algo == 'EA-vc' and rin.n_elim > 0:
            elimination_context = prepare_elimination_context(
                rin,
                ranking,
                selected_nat_data,
                cn_data=cn_data,
            )
            elim_parent_ids = list(elimination_context.elim_dnat_map)
            logger.info(
                f'eligible candidates for elimination: {len(elim_parent_ids)}'
            )
            selection_context_el = prepare_selection_context(
                rin,
                elim_parent_ids,
                fitness,
            )

        # ---------- substitution contexts for EA-vc
        substitution_context = None
        selection_context_sb = None
        if rin.algo == 'EA-vc' and rin.n_subs > 0:
            substitution_context = prepare_substitution_context(
                rin,
                mindist,
                ranking,
                selected_nat_data,
                cn_data=cn_data,
            )
            subs_parent_ids = list(substitution_context.subs_comb_map)
            logger.info(
                f'eligible candidates for substitution: {len(subs_parent_ids)}'
            )
            selection_context_sb = prepare_selection_context(
                rin,
                subs_parent_ids,
                fitness,
            )

        # ---------- offspring tasks
        tasks = create_offspring_tasks(
            conn,
            rin,
            next_generation,
        )

        # ---------- random IDs
        random_id_start = get_next_record_id(conn)
        random_ids = list(
            range(
                random_id_start,
                random_id_start + rin.n_rand,
            )
        )

        # ---------- population
        if len(tasks) + len(random_ids) != rin.n_pop:
            logger.error(
                'Number of offspring tasks and random IDs '
                f'does not match n_pop: '
                f'{len(tasks)} + {len(random_ids)} '
                f'!= {rin.n_pop}'
            )
            raise SystemExit(1)

        # ---------- offspring generator
        offspring_generator = prepare_offspring_generator(
            rin,
            parent_data,
            selection_context,
            selection_context_pm,
            crossover_context,
            permutation_context,
            strain_context,
            selection_context_ad,
            addition_context,
            selection_context_el,
            elimination_context,
            selection_context_sb,
            substitution_context,
        )

        # ---------- elite candidates
        elite_pool_struc = selected_struc_data
        elite_pool_fitness = fitness
        if rin.n_elite > 0:
            elite_pool_struc, elite_pool_fitness = prepare_elite_pool(
                conn,
                rin,
                current_generation,
            )

        # ---------- generate and save offspring
        elite_ids = generate_and_save_offspring_ea(
            conn,
            rin,
            tasks,
            offspring_generator,
            elite_pool_struc,
            elite_pool_fitness,
            rng,
        )

        # ---------- prepare random records
        generate_ids = prepare_random_records_ea(
            conn,
            random_ids,
            next_generation,
        )
        conn.commit()

        # ---------- generate random structures
        generate_structures_rs(
            rin,
            generate_ids,
            feasible_N=feasible_N,
            cn_data=cn_data,
        )

        try:
            # ---------- update random records
            update_generation(
                conn,
                random_ids,
                next_generation,
            )

            # ---------- save random origins
            for cid in random_ids:
                insert_ea_origin(
                    conn,
                    cid,
                    next_generation,
                    'random',
                )

            # ---------- save elite IDs
            save_ea_elite_ids(
                conn,
                next_generation,
                elite_ids,
            )

            # ---------- save generation information
            save_ea_generation_info(
                conn,
                rin,
                next_generation,
            )

            conn.commit()

        except Exception:
            # ---------- rollback
            conn.rollback()
            raise

    finally:
        # ---------- close database
        if conn is not None:
            conn.close()

    # ---------- time
    etime = timedelta(seconds=perf_counter() - time_start)
    logger.info(
        'Elapsed time for generation '
        f'{next_generation} structure generation: {etime}'
    )
