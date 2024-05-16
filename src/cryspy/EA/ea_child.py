from logging import getLogger

from .gen_struc_EA.crossover import Crossover
from .gen_struc_EA.ea_generation import EA_generation
from .gen_struc_EA.permutation import Permutation
from .gen_struc_EA.rotation import Rotation
from .gen_struc_EA.strain import Strain
from .gen_struc_EA.addition import Addition
from .gen_struc_EA.elimination import Elimination
from .gen_struc_EA.substitution import Substitution
from ..IO import pkl_data
from ..IO.out_results import out_nat_data
from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, get_mol_data, out_poscar, get_nat

# ---------- import later
#from ..RS.gen_struc_RS import gen_pyxtal
#from ..RS.gen_struc_RS import random_generation


logger = getLogger('cryspy')

def child_gen(rin, sp, init_struc_data, struc_mol_id=None, ea_vc_data=None):
    # ---------- instantiate EA_generation class
    eagen = EA_generation(sp=sp, id_start=rin.tot_struc)

    # ---------- set mindist
    logger.info('# mindist')
    mindist = set_mindist(rin.atype, rin.mindist, rin.mindist_factor, rin.struc_mode)
    if rin.struc_mode == 'mol_bs':
        mindist_dummy = set_mindist(rin.atype,
                                    rin.mindist_mol_bs,
                                    rin.mindist_mol_bs_factor,
                                    rin.struc_mode,
                                    dummy=True,
                                    mol_file=rin.mol_file,
                                    mpi_rank=0)

    # ------ Crossover
    if rin.n_crsov > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            co = Crossover(rin, mindist)
            eagen.gen_crossover(rin, co, None, molecular=False)
        else:
            co = Crossover(rin, mindist)
            eagen.gen_crossover(rin, co, struc_mol_id, molecular=True)

    # ------ Permutation
    if rin.n_perm > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            pm = Permutation(mindist)
            eagen.gen_permutation(rin, pm, None, molecular=False)
        else:
            pm = Permutation(mindist)
            eagen.gen_permutation(rin, pm, struc_mol_id, molecular=True)

    # ------ Strain
    if rin.n_strain > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            st = Strain(mindist)
            eagen.gen_strain(rin, st, None)
        else:
            st = Strain(mindist)
            eagen.gen_strain(rin, st, struc_mol_id, protect_mol_struc=True)

    # ------ EA-vc
    if rin.algo == 'EA-vc':
        nat_data, _, _ = ea_vc_data

        # -- Addition
        if rin.n_add > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                ad = Addition(mindist, rin.target)
                eagen.gen_addition(rin, ad, nat_data)
            else:
                raise SystemExit(1)

        # -- Elimination
        if rin.n_elim > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                el = Elimination(mindist, rin.target)
                eagen.gen_elimination(rin, el, nat_data)
            else:
                raise SystemExit(1)

        # -- Substitution
        if rin.n_subs > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                su = Substitution(mindist, rin.target)
                eagen.gen_substitution(rin, su, nat_data)
            else:
                raise SystemExit(1)

    # ------ Rotation
    if rin.struc_mode in ['mol', 'mol_bs']:
        if rin.n_rotation > 0:
            rot = Rotation(mindist)
            eagen.gen_rotation(rin, struc_mol_id, rot=rot)

    # ------ write init_POSCARS
    out_poscar(eagen.offspring, './data/init_POSCARS')

    # ------ update init_struc_data
    init_struc_data.update(eagen.offspring)
    if rin.struc_mode in ['mol', 'mol_bs']:
        struc_mol_id.update(eagen.offspring_mol_id)

    # ------ save EA-vc_data.pkl
    if rin.algo == 'EA-vc':
        nat_data, ratio_data, hdist_data = ea_vc_data
        for cid, struc in eagen.offspring.items():
            nat_data[cid], ratio_data[cid] = get_nat(struc, rin.atype)
        ea_vc_data = (nat_data, ratio_data, hdist_data)
        pkl_data.save_ea_vc_data(ea_vc_data)

    # ---------- random generation
    if rin.n_rand > 0:
        vc = True if rin.algo == 'EA-vc' else False
        # ------ pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            from ..RS.gen_struc_RS import gen_pyxtal
            # -- crystal
            if rin.struc_mode == 'crystal':
                tmp_struc_data = gen_pyxtal.gen_struc(rin,
                                                       nstruc=rin.n_rand,
                                                       mindist=mindist,
                                                       id_offset=eagen.cid,
                                                       vc=vc)
            # -- molecular crystal
            elif rin.struc_mode == 'mol':
                mol_data = get_mol_data(rin.mol_file)
                tmp_struc_data, tmp_mol_id = gen_pyxtal.gen_struc_mol(
                                                 rin,
                                                 nstruc=rin.n_rand,
                                                 mindist=mindist,
                                                 mol_data=mol_data,
                                                 id_offset=eagen.cid,
                                                 vc=False
                                             )
            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                logger.info('# -- mindist_mol_bs')
                mol_data = get_mol_data(rin.mol_file)
                tmp_struc_data, tmp_mol_id = gen_pyxtal.gen_struc_mol_break_sym(
                                                 rin,
                                                 nstruc=rin.n_rand,
                                                 mindist=mindist,
                                                 mindist_dummy=mindist_dummy,
                                                 mol_data=mol_data,
                                                 id_offset=eagen.cid
                                             )
        # ------ w/o pyxtal
        else:
            from ..RS.gen_struc_RS import random_generation
            if rin.spgnum == 0:
                init_struc_data = random_generation.gen_wo_spg(
                                      rin,
                                      nstruc=rin.n_rand,
                                      mindist=mindist,
                                      id_offset=eagen.cid,
                                      vc=vc
                                  )
            else:
                fwpath = check_fwpath(rin.fwpath)
                init_struc_data = random_generation.gen_with_find_wy(
                                    rin,
                                    nstruc=rin.n_rand,
                                    mindist=mindist,
                                    id_offset=eagen.cid,
                                    fwpath=fwpath,
                                    mpi_rank=0,
                                    vc=vc
                                )
        # ------ update
        init_struc_data.update(tmp_struc_data)
        if rin.struc_mode in ['mol', 'mol_bs']:
            struc_mol_id.update(tmp_mol_id)
        # ------ save EA-vc_data.pkl
        if rin.algo == 'EA-vc':
            # ea_vc_data is already loaded above
            for cid, struc in tmp_struc_data.items():
                nat_data[cid], ratio_data[cid] = get_nat(struc, rin.atype)
            ea_vc_data = (nat_data, ratio_data, hdist_data)
            pkl_data.save_ea_vc_data(ea_vc_data)
        # ------ write init_POSCARS
        out_poscar(tmp_struc_data, './data/init_POSCARS')

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)
    if rin.struc_mode in ['mol', 'mol_bs']:
        pkl_data.save_struc_mol_id(struc_mol_id)

    # ---------- out nat_data
    if rin.algo == 'EA-vc':
        out_nat_data(nat_data, rin.atype)

    # ----------return
    return init_struc_data, eagen, struc_mol_id


# not used in this version
# this is used in adj_comp.py
def check_vcnat(rin, child):
    from ..util.struc_util import get_nat
    nat, ratio = get_nat(child, rin.atype)

    if len(rin.atype) == 2:
        for i in range(len(rin.atype)):
            if (nat[i] >= rin.ll_nat[i]) and (nat[i]<= rin.ul_nat[i]):
                check_nat = True
            else:
                #check_nat = False
                return False

        return True
        
    elif len(rin.atype) == 3:
        SystemExit(1) #temporary
    else:
        SystemExit(1) #temporary