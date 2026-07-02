from logging import getLogger

from .gen_struc_EA.select_parents import SelectParents
from .gen_struc_EA.crossover import gen_crossover
from .gen_struc_EA.permutation import gen_permutation
from .gen_struc_EA.strain import gen_strain
from .gen_struc_EA.addition import gen_addition
from .gen_struc_EA.elimination import gen_elimination
from .gen_struc_EA.substitution import gen_substitution
#from .gen_struc_EA.rotation import gen_rotation
from ..IO import pkl_data
from ..IO.out_results import out_nat_data
from ..RS.rs_gen import gen_random
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
            mpi_rank=0,
        )
    # ---------- initialize
    children = {}
    parents = {}
    operation = {}
    pre_nstruc = len(init_struc_data)
    id_start = pre_nstruc
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

    # ---------- helper: build SelectParents with common selector settings
    def _build_sp(parent_ids):
        sp_local = SelectParents(parent_ids, rng)
        if rin.slct_func == 'TNM':
            sp_local.set_tournament(rin.t_size)
        else:
            # roulette needs fitness for IDs in parent_ids
            fit_local = {pid: fittest[pid] for pid in parent_ids}
            sp_local.set_roulette(fit_local, rin.a_rlt, rin.b_rlt, rin.fit_reverse)
            logger.debug(f'cumulative fitness in roulette: {sp_local.cum_fit}')
        return sp_local

    # ---------- parent selectors for each operation
    # ------ crossover and strain
    sp_co = _build_sp(ranking)
    logger.info(f'eligible candidates for crossover and strain: {len(sp_co.ranking)}')
    # ------ permutation
    if vc:
        sp_pm = _build_sp(perm_pid_list) if rin.n_perm > 0 else None
    else:
        sp_pm = _build_sp(ranking) if rin.n_perm > 0 else None
    if sp_pm is not None:
        logger.info(f'eligible candidates for permutation: {len(sp_pm.ranking)}')
    # ------ EA-vc only operation-specific parent pools
    sp_add = _build_sp(list(add_dnat_map.keys())) if (vc and rin.n_add > 0) else None
    sp_elim = _build_sp(list(elim_dnat_map.keys())) if (vc and rin.n_elim > 0) else None
    sp_subs = _build_sp(list(subs_comb_map.keys())) if (vc and rin.n_subs > 0) else None
    if sp_add is not None:
        logger.info(f'eligible candidates for addition: {len(sp_add.ranking)}')
    if sp_elim is not None:
        logger.info(f'eligible candidates for elimination: {len(sp_elim.ranking)}')
    if sp_subs is not None:
        logger.info(f'eligible candidates for substitution: {len(sp_subs.ranking)}')

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
            co_children, co_parents, co_operation = gen_crossover(
                rin.atype,
                rin.nat,
                mindist,
                struc_data,
                sp_co,
                rin.n_crsov,
                id_start,
                rin.symprec,
                rin.crs_lat,
                rin.nat_diff_tole,
                rin.maxcnt_ea,
                vc=vc,
                ll_nat=rin.ll_nat,
                ul_nat=rin.ul_nat,
                cn_comb=cn_comb_for_crossover,
                feasible_N=feasible_N,
                charge=rin.charge,
                cn_data=cn_data,
                min_comp=rin.min_comp,
                max_comp=rin.max_comp,
                rng=rng,
            )
        else:
            logger.error('Crossover is not implemented for mol or mol_bs')
        # ------ update
        children.update(co_children)
        parents.update(co_parents)
        operation.update(co_operation)
        id_start += rin.n_crsov

    # ---------- Permutation
    if rin.n_perm > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            pm_children, pm_parents, pm_operation = gen_permutation(
                rin.atype,
                mindist,
                struc_data,
                sp_pm,
                rin.n_perm,
                id_start,
                rin.symprec,
                rin.ntimes,
                rin.maxcnt_ea,
                rng=rng,
            )
        else:
            logger.error('Permutation is not implemented for mol or mol_bs')
        # ------ update
        children.update(pm_children)
        parents.update(pm_parents)
        operation.update(pm_operation)
        id_start += rin.n_perm

    # ---------- Strain
    if rin.n_strain > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            st_children, st_parents, st_operation = gen_strain(
                rin.atype,
                mindist,
                struc_data,
                sp_co,
                rin.n_strain,
                id_start,
                rin.symprec,
                rin.sigma_st,
                rin.maxcnt_ea,
                rng=rng,
            )
        else:
            logger.error('Strain is not implemented for mol or mol_bs')
        # ------ update
        children.update(st_children)
        parents.update(st_parents)
        operation.update(st_operation)
        id_start += rin.n_strain

    # ---------- EA-vc
    if rin.algo == 'EA-vc':

        # ------ Addition
        if rin.n_add > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                ad_children, ad_parents, ad_operation = gen_addition(
                    rin.atype,
                    mindist,
                    struc_data,
                    sp_add,
                    rin.n_add,
                    add_dnat_map,
                    id_start,
                    rin.symprec,
                    rin.maxcnt_ea,
                    rin.target,
                    rng=rng,
                )
            else:
                logger.error('Addition is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            children.update(ad_children)
            parents.update(ad_parents)
            operation.update(ad_operation)
            id_start += rin.n_add

        # ------ Elimination
        if rin.n_elim > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                el_children, el_parents, el_operation = gen_elimination(
                    rin.atype,
                    struc_data,
                    sp_elim,
                    rin.n_elim,
                    elim_dnat_map,
                    id_start,
                    rin.symprec,
                    rin.target,
                    rng=rng,
                )
            else:
                logger.error('Elimination is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            children.update(el_children)
            parents.update(el_parents)
            operation.update(el_operation)
            id_start += rin.n_elim

        # ------ Substitution
        if rin.n_subs > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                sb_children, sb_parents, sb_operation = gen_substitution(
                    rin.atype,
                    mindist,
                    struc_data,
                    sp_subs,
                    rin.n_subs,
                    subs_comb_map,
                    id_start,
                    rin.symprec,
                    rin.maxcnt_ea,
                    rin.target,
                    rng=rng,
                )
            else:
                logger.error('Substitution is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            children.update(sb_children)
            parents.update(sb_parents)
            operation.update(sb_operation)
            id_start += rin.n_subs

    # ---------- write init_POSCARS
    out_poscar(children, './data/init_POSCARS')

    # ---------- update init_struc_data
    init_struc_data.update(children)

    # ---------- save nat_data for EA-vc
    if rin.algo == 'EA-vc':
        for cid, struc in children.items():
            nat_data[cid] = get_nat(struc, rin.atype)
        pkl_data.save_nat_data(nat_data)

    # ---------- random generation
    if rin.n_rand > 0:
        logger.info('# ------ Random structure generation')
        tmp_struc_data = gen_random(
                                        rin=rin,
                                        nstruc=rin.n_rand,
                                        id_offset=id_start,
                                        comm=None,
                                        mpi_rank=0,
                                        mpi_size=1,
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

