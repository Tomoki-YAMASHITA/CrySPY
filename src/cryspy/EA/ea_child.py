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


logger = getLogger('cryspy')


def child_gen(
        rin,
        ranking,
        fittest,
        struc_data,
        init_struc_data,
        struc_mol_id=None,
        nat_data=None,
    ):

    # ---------- instantiate SelectParents class
    sp = SelectParents(ranking)    # after set_xxx, we can use sp.get_parents(n_parent)
    if rin.slct_func == 'TNM':
        sp.set_tournament(rin.t_size)
    else:
        sp.set_roulette(fittest, rin.a_rlt, rin.b_rlt, rin.fit_reverse)
        logger.debug(f'cumulative fitness in roulette: {sp.cum_fit}')

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
    #children_mol_id = {}
    parents = {}
    operation = {}
    pre_nstruc = len(init_struc_data)
    id_start = pre_nstruc
    vc = True if rin.algo == 'EA-vc' else False

    # ---------- Crossover
    if rin.n_crsov > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            co_children, co_parents, co_operation = gen_crossover(
                rin.atype,
                rin.nat,
                mindist,
                struc_data,
                sp,
                rin.n_crsov,
                id_start,
                rin.symprec,
                rin.crs_lat,
                rin.nat_diff_tole,
                rin.maxcnt_ea,
                vc,
                rin.ll_nat,
                rin.ul_nat,
                struc_mol_id=None,
                molecular=False,
            )
        else:
            logger.error('Crossover is not implemented for mol or mol_bs')
            # co = Crossover(rin, mindist)
            # eagen.gen_crossover(rin, co, struc_mol_id, molecular=True)
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
                sp,
                rin.n_perm,
                id_start,
                rin.symprec,
                rin.ntimes,
                rin.maxcnt_ea,
                struc_mol_id=None,
                molecular=False,
            )
        else:
            logger.error('Permutation is not implemented for mol or mol_bs')
            # pm = Permutation(mindist)
            # eagen.gen_permutation(rin, pm, struc_mol_id, molecular=True)
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
                sp,
                rin.n_strain,
                id_start,
                rin.symprec,
                rin.sigma_st,
                rin.maxcnt_ea,
                struc_mol_id=None,
                molecular=False,
                protect_mol_struc=True,
            )
        else:
            logger.error('Strain is not implemented for mol or mol_bs')
            # st = Strain(mindist)
            # eagen.gen_strain(rin, st, struc_mol_id, protect_mol_struc=True)
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
                    sp,
                    rin.n_add,
                    nat_data,
                    rin.ul_nat,
                    id_start,
                    rin.symprec,
                    rin.maxcnt_ea,
                    rin.target,
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
                    sp,
                    rin.n_elim,
                    nat_data,
                    rin.ll_nat,
                    id_start,
                    rin.symprec,
                    rin.target,
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
                    sp,
                    rin.n_subs,
                    nat_data,
                    rin.ll_nat,
                    rin.ul_nat,
                    id_start,
                    rin.symprec,
                    rin.maxcnt_ea,
                    rin.target,
                )
            else:
                logger.error('Substitution is not implemented for mol or mol_bs')
                raise SystemExit(1)
            # ------ update
            children.update(sb_children)
            parents.update(sb_parents)
            operation.update(sb_operation)
            id_start += rin.n_subs

    # ---------- Rotation
    # if rin.struc_mode in ['mol', 'mol_bs']:
    #     if rin.n_rotation > 0:
    #         rot = Rotation(mindist)
    #         eagen.gen_rotation(rin, struc_mol_id, rot=rot)

    # ---------- write init_POSCARS
    out_poscar(children, './data/init_POSCARS')

    # ---------- update init_struc_data
    init_struc_data.update(children)
    # if rin.struc_mode in ['mol', 'mol_bs']:
    #     struc_mol_id.update(children_mol_id)

    # ---------- save EA-vc_data.pkl
    if rin.algo == 'EA-vc':
        for cid, struc in children.items():
            nat_data[cid] = get_nat(struc, rin.atype)
        pkl_data.save_nat_data(nat_data)

    # ---------- random generation
    if rin.n_rand > 0:
        tmp_struc_data, tmp_mol_id = gen_random(
            rin,
            rin.n_rand,
            id_start,
            comm=None,
            mpi_rank=0,
            mpi_size=1,
        )
        # ------ update
        init_struc_data.update(tmp_struc_data)
        # if rin.struc_mode in ['mol', 'mol_bs']:
        #     struc_mol_id.update(tmp_mol_id)
        # ------ save EA-vc_data.pkl
        if rin.algo == 'EA-vc':
            for cid, struc in tmp_struc_data.items():
                nat_data[cid] = get_nat(struc, rin.atype)
            pkl_data.save_nat_data(nat_data)
        # ------ write init_POSCARS
        out_poscar(tmp_struc_data, './data/init_POSCARS')

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)
    # if rin.struc_mode in ['mol', 'mol_bs']:
    #     pkl_data.save_struc_mol_id(struc_mol_id)

    # ---------- out nat_data
    if rin.algo == 'EA-vc':
        out_nat_data(nat_data, rin.atype)

    # ----------return
    return init_struc_data, parents, operation
    #return init_struc_data, parents, operation, struc_mol_id
