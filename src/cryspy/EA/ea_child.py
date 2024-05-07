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
from ..IO import read_input as rin
from ..IO.out_results import out_nat_data
from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, out_poscar, get_nat

# ---------- import later
#from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
#from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen


logger = getLogger('cryspy')

def child_gen(sp, init_struc_data, struc_mol_id=None, ea_vc_data=None):
    # ---------- instantiate EA_generation class
    eagen = EA_generation(sp=sp, id_start=rin.tot_struc,
                          init_pos_path='./data/init_POSCARS')

    # ---------- set mindist
    logger.info('# mindist')
    mindist = set_mindist(rin.mindist, rin.mindist_factor)
    if rin.struc_mode == 'mol_bs':
        mindist_dummy = set_mindist(rin.mindist_mol_bs, rin.mindist_mol_bs_factor)

    # ------ Crossover
    if rin.n_crsov > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            co = Crossover(mindist)
            eagen.gen_crossover(co, None, molecular=False)
        else:
            co = Crossover(mindist)
            eagen.gen_crossover(co, struc_mol_id, molecular=True)

    # ------ Permutation
    if rin.n_perm > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            pm = Permutation(mindist)
            eagen.gen_permutation(pm, None, molecular=False)
        else:
            pm = Permutation(mindist)
            eagen.gen_permutation(pm, struc_mol_id, molecular=True)

    # ------ Strain
    if rin.n_strain > 0:
        if rin.struc_mode not in ['mol', 'mol_bs']:
            st = Strain(mindist)
            eagen.gen_strain(st, None)
        else:
            st = Strain(mindist)
            eagen.gen_strain(st, struc_mol_id, protect_mol_struc=True)

    # ------ EA-vc
    if rin.algo == 'EA-vc':
        nat_data, _, _ = ea_vc_data

        # -- Addition
        if rin.n_add > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                ad = Addition(mindist, rin.target)
                eagen.gen_addition(ad, nat_data)
            else:
                raise SystemExit(1)

        # -- Elimination
        if rin.n_elim > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                el = Elimination(mindist, rin.target)
                eagen.gen_elimination(el, nat_data)
            else:
                raise SystemExit(1)

        # -- Substitution
        if rin.n_subs > 0:
            if rin.struc_mode not in ['mol', 'mol_bs']:
                su = Substitution(mindist, rin.target)
                eagen.gen_substitution(su, nat_data)
            else:
                raise SystemExit(1)

    # ------ Rotation
    if rin.struc_mode in ['mol', 'mol_bs']:
        if rin.n_rotation > 0:
            rot = Rotation(mindist)
            eagen.gen_rotation(struc_mol_id, rot=rot)

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
        # ------ pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
            rsgx = Rnd_struc_gen_pyxtal(mindist=mindist)
            # -- crystal
            if rin.struc_mode == 'crystal':
                if not rin.algo == 'EA-vc':
                    rsgx.gen_struc(nstruc=rin.n_rand, id_offset=eagen.cid)
                else:    # vc
                    rsgx.gen_struc(nstruc=rin.n_rand, id_offset=eagen.cid, vc=True)
            # -- molecular crystal
            elif rin.struc_mode == 'mol':
                rsgx.set_mol()
                rsgx.gen_struc_mol(nstruc=rin.n_rand, id_offset=eagen.cid)
            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                logger.info('# -- mindist_mol_bs')
                mindist_dummy = set_mindist(rin.mindist_mol_bs, rin.mindist_mol_bs_factor, dummy=True)
                rsgx.set_mol()
                rsgx.gen_struc_mol_break_sym(nstruc=rin.n_rand,
                                             mindist_dummy=mindist_dummy,
                                             id_offset=eagen.cid)
            # -- update
            init_struc_data.update(rsgx.init_struc_data)
            if rin.struc_mode in ['mol', 'mol_bs']:
                struc_mol_id.update(rsgx.struc_mol_id)
            # -- save EA-vc_data.pkl
            if rin.algo == 'EA-vc':
                # ea_vc_data is already loaded above
                for cid, struc in rsgx.init_struc_data.items():
                    nat_data[cid], ratio_data[cid] = get_nat(struc, rin.atype)
                ea_vc_data = (nat_data, ratio_data, hdist_data)
                pkl_data.save_ea_vc_data(ea_vc_data)
            # -- init_POSCARS
            for cid, struc in rsgx.init_struc_data.items():
                out_poscar(struc, cid, './data/init_POSCARS')
        # ------ Rnd_struc_gen
        else:
            from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen
            rsg = Rnd_struc_gen(mindist=mindist)
            if rin.spgnum == 0:
                if not rin.algo == 'EA-vc':
                    rsg.gen_wo_spg(nstruc=rin.n_rand, id_offset=eagen.cid)
                else:    # vc
                    rsg.gen_wo_spg(nstruc=rin.n_rand, id_offset=eagen.cid, vc=True)
            else:
                fwpath = check_fwpath(rin.fwpath)
                if not rin.algo == 'EA-vc':
                    rsg.gen_with_find_wy(nstruc=rin.n_rand,
                                         id_offset=eagen.cid,
                                         fwpath=fwpath)
                else:
                    rsg.gen_with_find_wy(nstruc=rin.n_rand,
                                         id_offset=eagen.cid,
                                         fwpath=fwpath, vc=True)
            # -- update
            init_struc_data.update(rsg.init_struc_data)
            # -- save EA-vc_data.pkl
            if rin.algo == 'EA-vc':
                # ea_vc_data is already loaded above
                for cid, struc in rsg.init_struc_data.items():
                    nat_data[cid], ratio_data[cid] = get_nat(struc, rin.atype)
                ea_vc_data = (nat_data, ratio_data, hdist_data)
                pkl_data.save_ea_vc_data(ea_vc_data)
            # -- init_POSCARS
            for cid, struc in rsg.init_struc_data.items():
                out_poscar(struc, cid, './data/init_POSCARS')

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
def check_vcnat(child):
    from ..util.struc_util import get_nat
    nat, ratio = get_nat(child,rin.atype)

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