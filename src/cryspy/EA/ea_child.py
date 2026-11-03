from logging import getLogger

from .ea_offspring import (
    ParentData,
    TournamentSelectionContext,
    build_roulette_selection_context,
)
from .gen_struc_EA.crossover import (
    CrossoverContext,
    CrossoverOffspringGenerator,
    gen_crossover_batch,
)
from .gen_struc_EA.permutation import (
    PermutationContext,
    PermutationOffspringGenerator,
    gen_permutation_batch,
)
from .gen_struc_EA.strain import (
    StrainContext,
    StrainOffspringGenerator,
    gen_strain_batch,
)
from .gen_struc_EA.addition import (
    AdditionContext,
    AdditionOffspringGenerator,
    gen_addition_batch,
)
from .gen_struc_EA.elimination import (
    EliminationContext,
    EliminationOffspringGenerator,
    gen_elimination_batch,
)
from .gen_struc_EA.substitution import (
    SubstitutionContext,
    SubstitutionOffspringGenerator,
    gen_substitution_batch,
)
#from .gen_struc_EA.rotation import gen_rotation
from ..IO import pkl_data
from ..IO.out_results import out_nat_data
from ..RS.rs_gen import gen_random_batch
from ..util.struc_util import set_mindist, out_poscar, get_nat
#from ..util.struc_util import get_mol_data

# ---------- import later
#from .gen_struc_EA.addition import get_add_dnat_comb
#from .gen_struc_EA.elimination import get_elim_dnat_comb
#from .gen_struc_EA.substitution import get_subs_comb
#from ..util.struc_util import get_feasible_composition, precompute_feasible_N
#from ..util.charge_neutral import filter_cn_comb_comp


logger = getLogger('cryspy')


def child_gen(
        rin,
        ranking,
        fittest,
        struc_data,
        init_struc_data,
        nat_data=None,
        rng=None,
    ):

    # ---------- set mindist
    logger.info('# -- mindist')
    mindist = set_mindist(rin.atype, rin.mindist, rin.mindist_factor, rin.struc_mode)
    if rin.struc_mode == 'mol_bs':
        mindist_dummy = set_mindist(
            rin.atype,
            rin.mindist_mol_bs,
            rin.mindist_mol_bs_factor,
            rin.struc_mode,
            dummy=True,
            mol_file=rin.mol_file,
        )
    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}
    id_start = len(init_struc_data) + 1
    vc = True if rin.algo == 'EA-vc' else False
    use_comp_window = (rin.min_comp is not None or rin.max_comp is not None)
    use_charge_neutral = (rin.charge is not None)
    # ------ vc: charge neutrality
    # cn_data keeps the charge-neutral metadata loaded from cn_comb_data.pkl.
    # cn_comb is used only in enumerate mode; sample mode keeps cn_comb as None.
    cn_data = None
    cn_comb = None
    if vc and use_charge_neutral:
        cn_data = pkl_data.load_cn_comb_data()
        logger.info(f'Charge-neutral mode for child generation: {cn_data["mode"]}')

        # -- enumerate mode: use precomputed charge-neutral nat combinations
        if cn_data.get('mode') == 'enumerate':
            cn_comb = cn_data['cn_comb']

        # -- sample mode: cn_comb is not available; use cn_data and charge downstream
        elif cn_data.get('mode') == 'sample':
            cn_comb = None
        # -- unexpected
        else:
            logger.error('cn_data["mode"] must be enumerate or sample.')
            raise SystemExit(1)


    # ------ vc: pre-check parent feasibility and precompute per-parent operation data
    perm_pid_list = None
    if vc:
        perm_pid_list, add_dnat_map, elim_dnat_map, subs_comb_map = _check_parent_feasibility(
            rin, ranking, nat_data, cn_comb=cn_comb
        )

    # ---------- parent data
    parent_data = ParentData(
        structures=struc_data,
    )

    # ---------- log eligible candidates for each operation
    # ------ crossover and strain
    logger.info(f'eligible candidates for crossover and strain: {len(ranking)}')
    # ------ permutation
    if rin.n_perm > 0:
        perm_parent_ids = perm_pid_list if vc else ranking
        logger.info(f'eligible candidates for permutation: {len(perm_parent_ids)}')
    # ---------- selection contexts
    if rin.slct_func == 'TNM':
        selection_context_co = TournamentSelectionContext(
            parent_ids=tuple(ranking),
            tournament_size=rin.t_size,
        )
        selection_context_pm = (
            TournamentSelectionContext(
                parent_ids=tuple(perm_pid_list if vc else ranking),
                tournament_size=rin.t_size,
            )
            if rin.n_perm > 0
            else None
        )
    else:
        selection_context_co = build_roulette_selection_context(
            parent_ids=tuple(ranking),
            fitness=fittest,
            a_rlt=rin.a_rlt,
            b_rlt=rin.b_rlt,
            fit_reverse=rin.fit_reverse,
        )
        selection_context_pm = (
            build_roulette_selection_context(
                parent_ids=tuple(perm_pid_list if vc else ranking),
                fitness=fittest,
                a_rlt=rin.a_rlt,
                b_rlt=rin.b_rlt,
                fit_reverse=rin.fit_reverse,
            )
            if rin.n_perm > 0
            else None
        )

    # ---------- selection contexts for EA-vc
    # ------ selection context for addition
    selection_context_ad = None
    if vc and rin.n_add > 0:
        add_parent_ids = tuple(add_dnat_map)
        logger.info(f'eligible candidates for addition: {len(add_parent_ids)}')
        if rin.slct_func == 'TNM':
            selection_context_ad = TournamentSelectionContext(
                parent_ids=add_parent_ids,
                tournament_size=rin.t_size,
            )
        else:
            selection_context_ad = build_roulette_selection_context(
                parent_ids=add_parent_ids,
                fitness=fittest,
                a_rlt=rin.a_rlt,
                b_rlt=rin.b_rlt,
                fit_reverse=rin.fit_reverse,
            )

    # ------ selection context for elimination
    selection_context_el = None
    if vc and rin.n_elim > 0:
        elim_parent_ids = tuple(elim_dnat_map)
        logger.info(f'eligible candidates for elimination: {len(elim_parent_ids)}')
        if rin.slct_func == 'TNM':
            selection_context_el = TournamentSelectionContext(
                parent_ids=elim_parent_ids,
                tournament_size=rin.t_size,
            )
        else:
            selection_context_el = build_roulette_selection_context(
                parent_ids=elim_parent_ids,
                fitness=fittest,
                a_rlt=rin.a_rlt,
                b_rlt=rin.b_rlt,
                fit_reverse=rin.fit_reverse,
            )

    # ------ selection context for substitution
    selection_context_sb = None
    if vc and rin.n_subs > 0:
        subs_parent_ids = tuple(subs_comb_map)
        logger.info(f'eligible candidates for substitution: {len(subs_parent_ids)}')
        if rin.slct_func == 'TNM':
            selection_context_sb = TournamentSelectionContext(
                parent_ids=subs_parent_ids,
                tournament_size=rin.t_size,
            )
        else:
            selection_context_sb = build_roulette_selection_context(
                parent_ids=subs_parent_ids,
                fitness=fittest,
                a_rlt=rin.a_rlt,
                b_rlt=rin.b_rlt,
                fit_reverse=rin.fit_reverse,
            )

    # ---------- vc: prepare composition data for crossover and random generation
    feasible_comp = None
    feasible_N = None
    cn_comb_for_crossover = cn_comb

    if vc and use_comp_window:
        from ..util.struc_util import get_feasible_composition, precompute_feasible_N
        feasible_comp = get_feasible_composition(rin.min_comp, rin.max_comp)
        feasible_N = precompute_feasible_N(rin.ll_nat, rin.ul_nat, feasible_comp)
        logger.info(f'Feasible total atom counts for composition constraints: {len(feasible_N)}')

        # ------ enumerate mode: filter charge-neutral combinations by composition window
        if use_charge_neutral and cn_data['mode'] == 'enumerate':
            from ..util.charge_neutral import filter_cn_comb_comp
            cn_comb_for_crossover = filter_cn_comb_comp(
                cn_comb,
                min_comp=rin.min_comp,
                max_comp=rin.max_comp,
            )
            if len(cn_comb_for_crossover) == 0:
                logger.error('No charge-neutral compositions satisfy min/max_comp.')
                raise SystemExit(1)

    # ---------- Crossover
    if rin.n_crsov > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            co_generator = CrossoverOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context_co,
                context=CrossoverContext(
                    atype=rin.atype,
                    nat=rin.nat,
                    mindist=mindist,
                    crs_lat=rin.crs_lat,
                    nat_diff_tole=rin.nat_diff_tole,
                    maxcnt_ea=rin.maxcnt_ea,
                    vc=vc,
                    ll_nat=rin.ll_nat,
                    ul_nat=rin.ul_nat,
                    cn_comb=cn_comb_for_crossover,
                    feasible_N=feasible_N,
                    charge=rin.charge,
                    cn_data=cn_data,
                    min_comp=rin.min_comp,
                    max_comp=rin.max_comp,
                ),
                max_parent_attempts=rin.max_parent_attempts,
            )
            co_offspring_data, co_parents, co_operation = gen_crossover_batch(
                generator=co_generator,
                n_crsov=rin.n_crsov,
                id_start=id_start,
                symprec=rin.symprec,
                atype=rin.atype,
                base_seed=rin.seed,
                rng=rng,
            )
        else:
            logger.error('Crossover is not implemented for mol or mol_bs')
        # ------ update
        offspring_data.update(co_offspring_data)
        parents.update(co_parents)
        operation.update(co_operation)
        id_start += rin.n_crsov

    # ---------- Permutation
    if rin.n_perm > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            pm_generator = PermutationOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context_pm,
                context=PermutationContext(
                    atype=rin.atype,
                    mindist=mindist,
                    ntimes=rin.ntimes,
                    maxcnt_ea=rin.maxcnt_ea,
                ),
                max_parent_attempts=rin.max_parent_attempts,
            )
            pm_offspring_data, pm_parents, pm_operation = gen_permutation_batch(
                generator=pm_generator,
                n_perm=rin.n_perm,
                id_start=id_start,
                symprec=rin.symprec,
                atype=rin.atype,
                base_seed=rin.seed,
                rng=rng,
            )
        else:
            logger.error('Permutation is not implemented for mol or mol_bs')
        # ------ update
        offspring_data.update(pm_offspring_data)
        parents.update(pm_parents)
        operation.update(pm_operation)
        id_start += rin.n_perm

    # ---------- Strain
    if rin.n_strain > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            st_generator = StrainOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context_co,
                context=StrainContext(
                    atype=rin.atype,
                    mindist=mindist,
                    sigma_st=rin.sigma_st,
                    maxcnt_ea=rin.maxcnt_ea,
                ),
                max_parent_attempts=rin.max_parent_attempts,
            )
            st_offspring_data, st_parents, st_operation = gen_strain_batch(
                generator=st_generator,
                n_strain=rin.n_strain,
                id_start=id_start,
                symprec=rin.symprec,
                atype=rin.atype,
                base_seed=rin.seed,
                rng=rng,
            )
        else:
            logger.error('Strain is not implemented for mol or mol_bs')
        # ------ update
        offspring_data.update(st_offspring_data)
        parents.update(st_parents)
        operation.update(st_operation)
        id_start += rin.n_strain

    # ---------- EA-vc
    if rin.algo == 'EA-vc':

        # ------ Addition
        if rin.n_add > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                addition_context = AdditionContext(
                    atype=rin.atype,
                    mindist=mindist,
                    add_dnat_map=add_dnat_map,
                    maxcnt_ea=rin.maxcnt_ea,
                )
                ad_generator = AdditionOffspringGenerator(
                    parent_data=parent_data,
                    selection_context=selection_context_ad,
                    context=addition_context,
                    max_parent_attempts=rin.max_parent_attempts,
                )
                ad_children, ad_parents, ad_operation = gen_addition_batch(
                    generator=ad_generator,
                    n_add=rin.n_add,
                    id_start=id_start,
                    symprec=rin.symprec,
                    atype=rin.atype,
                    base_seed=rin.seed,
                    rng=rng,
                )
            else:
                logger.error('Addition is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            offspring_data.update(ad_children)
            parents.update(ad_parents)
            operation.update(ad_operation)
            id_start += rin.n_add

        # ------ Elimination
        if rin.n_elim > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                elimination_context = EliminationContext(
                    atype=rin.atype,
                    elim_dnat_map=elim_dnat_map,
                )
                el_generator = EliminationOffspringGenerator(
                    parent_data=parent_data,
                    selection_context=selection_context_el,
                    context=elimination_context,
                    max_parent_attempts=rin.max_parent_attempts,
                )
                el_children, el_parents, el_operation = gen_elimination_batch(
                    generator=el_generator,
                    n_elim=rin.n_elim,
                    id_start=id_start,
                    symprec=rin.symprec,
                    atype=rin.atype,
                    base_seed=rin.seed,
                    rng=rng,
                )
            else:
                logger.error('Elimination is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            offspring_data.update(el_children)
            parents.update(el_parents)
            operation.update(el_operation)
            id_start += rin.n_elim

        # ------ Substitution
        if rin.n_subs > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                substitution_context = SubstitutionContext(
                    atype=rin.atype,
                    mindist=mindist,
                    subs_comb_map=subs_comb_map,
                    maxcnt_ea=rin.maxcnt_ea,
                )
                sb_generator = SubstitutionOffspringGenerator(
                    parent_data=parent_data,
                    selection_context=selection_context_sb,
                    context=substitution_context,
                    max_parent_attempts=rin.max_parent_attempts,
                )
                sb_children, sb_parents, sb_operation = gen_substitution_batch(
                    generator=sb_generator,
                    n_subs=rin.n_subs,
                    id_start=id_start,
                    symprec=rin.symprec,
                    atype=rin.atype,
                    base_seed=rin.seed,
                    rng=rng,
                )
            else:
                logger.error('Substitution is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            offspring_data.update(sb_children)
            parents.update(sb_parents)
            operation.update(sb_operation)
            id_start += rin.n_subs

    # ---------- write init_POSCARS
    out_poscar(offspring_data, './data/init_POSCARS')

    # ---------- update init_struc_data
    init_struc_data.update(offspring_data)

    # ---------- save nat_data for EA-vc
    if rin.algo == 'EA-vc':
        for cid, struc in offspring_data.items():
            nat_data[cid] = get_nat(struc, rin.atype)
        pkl_data.save_nat_data(nat_data)

    # ---------- random generation
    if rin.n_rand > 0:
        logger.info('# ------ Random structure generation')
        tmp_struc_data = gen_random_batch(
                                        rin=rin,
                                        nstruc=rin.n_rand,
                                        id_offset=id_start,
                                        cn_data=cn_data,
                                        feasible_N=feasible_N,
                                        rng=rng,
                                    )
        # ------ update
        init_struc_data.update(tmp_struc_data)
        # ------ save nat_data for EA-vc
        if rin.algo == 'EA-vc':
            for cid, struc in tmp_struc_data.items():
                nat_data[cid] = get_nat(struc, rin.atype)
            pkl_data.save_nat_data(nat_data)
        # ------ write init_POSCARS
        out_poscar(tmp_struc_data, './data/init_POSCARS')

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)

    # ---------- out nat_data
    if rin.algo == 'EA-vc':
        out_nat_data(nat_data, rin.atype)

    # ----------return
    return init_struc_data, parents, operation


def _check_parent_feasibility(rin, ranking, nat_data, cn_comb=None, tol=1e-12):
    """
    Validate parent feasibility and precompute per-parent operation data.

    Returns
    -------
    perm_pid_list : list[int]
        Parent IDs that can be used for permutation (>=2 non-zero elements).
    add_dnat_map : dict[int, list[tuple[int, ...]]]
        Per-parent addition candidates: `add_dnat_map[pid] = dnat_comb`.
        (`dnat_comb` is the list of feasible delta-nat tuples for that parent.)
    elim_dnat_map : dict[int, list[tuple[int, ...]]]
        Per-parent elimination candidates: `elim_dnat_map[pid] = dnat_comb`.
    subs_comb_map : dict[int, list[list[tuple[str, str]]]]
        Per-parent substitution candidates: `subs_comb_map[pid] = subs_comb`.
    """
    def abort(msg):
        logger.error(msg)
        raise SystemExit(1)

    def in_comp_range(nat):
        if rin.min_comp is None or rin.max_comp is None:
            return True
        ntot = sum(nat)
        if ntot <= 0:
            return False
        return all(
            (cmin - tol) <= (n_i / ntot) <= (cmax + tol)
            for n_i, cmin, cmax in zip(nat, rin.min_comp, rin.max_comp)
        )

    def apply_subs_to_nat(parent_nat, subs_list, atype):
        idx = {a: i for i, a in enumerate(atype)}
        new_nat = list(parent_nat)
        for from_elem, to_elem in subs_list:
            i = idx[from_elem]
            j = idx[to_elem]
            new_nat[i] -= 1
            new_nat[j] += 1
            if new_nat[i] < 0:
                return None
        return tuple(new_nat)

    from .gen_struc_EA.addition import get_add_dnat_comb
    from .gen_struc_EA.elimination import get_elim_dnat_comb
    from .gen_struc_EA.substitution import get_subs_comb

    # ---------- valid parent IDs from ranking
    ranked_pids = []
    for cid in ranking:
        nat = nat_data.get(cid)
        if nat is not None and sum(nat) > 0:
            ranked_pids.append(cid)

    # ---------- crossover feasibility
    if rin.n_crsov > 0 and len(ranked_pids) < 2:
        abort(
            f'Crossover requires at least 2 parent candidates, but got {len(ranked_pids)}.'
        )

    # ---------- permutation feasibility + list
    perm_pid_list = [pid for pid in ranked_pids if sum(1 for n in nat_data[pid] if n > 0) >= 2]
    if rin.n_perm > 0 and len(perm_pid_list) == 0:
        abort('Permutation requires at least one parent containing >=2 non-zero elements.')

    # ---------- default maps (for non EA-vc or disabled ops)
    add_dnat_map = None
    elim_dnat_map = None
    subs_comb_map = None

    # ---------- EA-vc only: precompute add/elim/subs feasibility data
    if rin.algo == 'EA-vc':
        # ------ addition
        if rin.n_add > 0:
            add_dnat_map = {}
            for pid in ranked_pids:
                parent_nat = tuple(nat_data[pid])
                dnat_list = get_add_dnat_comb(
                    rin.ul_nat,
                    rin.add_max,
                    parent_nat,
                    cn_comb=cn_comb,
                    charge=rin.charge,
                )
                dnat_ok = []
                for dnat in dnat_list:
                    new_nat = tuple(a + b for a, b in zip(parent_nat, dnat))
                    if in_comp_range(new_nat):
                        dnat_ok.append(tuple(dnat))
                if dnat_ok:
                    add_dnat_map[pid] = dnat_ok
            if len(add_dnat_map) == 0:
                abort('Addition has no feasible parent under current nat/composition constraints.')

        # ------ elimination
        if rin.n_elim > 0:
            elim_dnat_map = {}
            for pid in ranked_pids:
                parent_nat = tuple(nat_data[pid])
                dnat_list = get_elim_dnat_comb(
                    rin.ll_nat,
                    rin.elim_max,
                    parent_nat,
                    cn_comb=cn_comb,
                    charge=rin.charge,
                )
                dnat_ok = []
                for dnat in dnat_list:
                    new_nat = tuple(a - b for a, b in zip(parent_nat, dnat))
                    if sum(new_nat) > 0 and in_comp_range(new_nat):
                        dnat_ok.append(tuple(dnat))
                if dnat_ok:
                    elim_dnat_map[pid] = dnat_ok
            if len(elim_dnat_map) == 0:
                abort('Elimination has no feasible parent under current nat/composition constraints.')

        # ------ substitution
        if rin.n_subs > 0:
            subs_comb_map = {}
            cn_set = {tuple(row) for row in cn_comb} if cn_comb is not None else None
            for pid in ranked_pids:
                parent_nat = tuple(nat_data[pid])
                subs_patterns = get_subs_comb(
                    rin.atype,
                    rin.ll_nat,
                    rin.ul_nat,
                    rin.subs_max,
                    parent_nat,
                    cn_set=cn_set,
                    charge=rin.charge,
                )
                subs_ok = []
                for subs_list in subs_patterns:
                    new_nat = apply_subs_to_nat(parent_nat, subs_list, rin.atype)
                    if new_nat is not None and in_comp_range(new_nat):
                        subs_ok.append(subs_list)
                if subs_ok:
                    subs_comb_map[pid] = subs_ok
            if len(subs_comb_map) == 0:
                abort('Substitution has no feasible parent under current nat/composition constraints.')

    # ---------- return
    return perm_pid_list, add_dnat_map, elim_dnat_map, subs_comb_map

