from logging import getLogger

from .gen_struc_EA.crossover import Crossover
from .gen_struc_EA.ea_generation import EA_generation
from .gen_struc_EA.permutation import Permutation
from .gen_struc_EA.rotation import Rotation
from .gen_struc_EA.strain import Strain
from ..IO import pkl_data
from ..IO import read_input as rin
from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, out_poscar

# ---------- import later
#from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
#from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen


logger = getLogger('cryspy')

def child_gen(sp, init_struc_data, struc_mol_id=None):
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

    # ------ Rotation
    if rin.struc_mode in ['mol', 'mol_bs']:
        if rin.n_rotation > 0:
            rot = Rotation(mindist)
            eagen.gen_rotation(struc_mol_id, rot=rot)

    # ------ update init_struc_data
    init_struc_data.update(eagen.offspring)
    if rin.struc_mode in ['mol', 'mol_bs']:
        struc_mol_id.update(eagen.offspring_mol_id)

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
            # -- init_POSCARS
            for cid, struc in rsg.init_struc_data.items():
                out_poscar(struc, cid, './data/init_POSCARS')

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)
    if rin.struc_mode in ['mol', 'mol_bs']:
        pkl_data.save_struc_mol_id(struc_mol_id)

    # ----------return
    if rin.struc_mode not in ['mol', 'mol_bs']:
        return init_struc_data, eagen
    else:
        return init_struc_data, eagen, struc_mol_id
