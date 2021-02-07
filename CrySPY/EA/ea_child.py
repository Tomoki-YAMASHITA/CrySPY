from .. import utility
from ..gen_struc.EA.ea_generation import EA_generation
from ..gen_struc.EA.crossover import Crossover
from ..gen_struc.EA.permutation import Permutation
from ..gen_struc.EA.strain import Strain
from ..gen_struc.random.random_generation import Rnd_struc_gen
from ..gen_struc.random.gen_pyxtal import Rnd_struc_gen_pyxtal
from ..IO import pkl_data
from ..IO import read_input as rin


def child_gen(sp, init_struc_data):
    # ---------- instantiate EA_generation class
    eagen = EA_generation(sp=sp, symprec=rin.symprec, id_start=rin.tot_struc,
                          init_pos_path='./data/init_POSCARS')

    # ------ instantiate Crossover class
    if rin.n_crsov > 0:
        co = Crossover(rin.atype, rin.nat, rin.mindist_ea,
                       rin.crs_lat, rin.nat_diff_tole, rin.maxcnt_ea)
        eagen.gen_crossover(rin.n_crsov, co=co)    # crossover
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by crossover\n'.format(rin.n_crsov))

    # ------ instantiate Permutation class
    if rin.n_perm > 0:
        pm = Permutation(rin.atype, rin.mindist_ea, rin.ntimes, rin.maxcnt_ea)
        eagen.gen_permutation(rin.n_perm, pm=pm)    # permutation
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by permutation\n'.format(rin.n_perm))

    # ------ instantiate Strain class
    if rin.n_strain > 0:
        st = Strain(rin.atype, rin.mindist_ea, rin.sigma_st, rin.maxcnt_ea)
        eagen.gen_strain(rin.n_strain, st=st)    # strain
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by strain\n'.format(rin.n_strain))

    # ------ update init_struc_data
    init_struc_data.update(eagen.offspring)

    # ---------- random generation
    if rin.n_rand > 0:
        # ------ pyxtal
        if not (rin.spgnum == 0 or rin.use_find_wy):
            rsgx = Rnd_struc_gen_pyxtal(natot=rin.natot, atype=rin.atype,
                                        nat=rin.nat, vol_factor=rin.vol_factor,
                                        vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                                        mindist=rin.mindist,
                                        spgnum=rin.spgnum, symprec=rin.symprec)
            # -- crystal
            if rin.struc_mode == 'crystal':
                rsgx.gen_struc(nstruc=rin.n_rand, id_offset=eagen.cid,
                               init_pos_path='./data/init_POSCARS')
            # -- molecular crystal
            elif rin.struc_mode == 'mol':
                rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
                rsgx.gen_struc_mol(nstruc=rin.n_rand, id_offset=eagen.cid,
                                   init_pos_path='./data/init_POSCARS',
                                   timeout_mol=rin.timeout_mol)
            # ------ molecular crystal breaking symmetry
            elif rin.struc_mode == 'mol_bs':
                rsgx.set_mol(mol_file=rin.mol_file, nmol=rin.nmol)
                rsgx.gen_struc_mol_break_sym(nstruc=rin.n_rand,
                                             rot_mol=rin.rot_mol,
                                             id_offset=eagen.cid,
                                             init_pos_path='./data/init_POSCARS')
            # -- update
            init_struc_data.update(rsgx.init_struc_data)
        # ------ Rnd_struc_gen
        else:
            rsg = Rnd_struc_gen(natot=rin.natot, atype=rin.atype, nat=rin.nat,
                                minlen=rin.minlen, maxlen=rin.maxlen,
                                dangle=rin.dangle, mindist=rin.mindist,
                                vol_mu=rin.vol_mu, vol_sigma=rin.vol_sigma,
                                maxcnt=rin.maxcnt, symprec=rin.symprec)
            if rin.spgnum == 0:
                rsg.gen_wo_spg(nstruc=rin.n_rand, id_offset=eagen.cid,
                               init_pos_path='./data/init_POSCARS')
                init_struc_data.update(rsg.init_struc_data)
            else:
                fwpath = utility.check_fwpath()
                rsg.gen_with_find_wy(nstruc=rin.n_rand, spgnum=rin.spgnum,
                                     id_offset=eagen.cid,
                                     init_pos_path='./data/init_POSCARS',
                                     fwpath=fwpath)
                init_struc_data.update(rsg.init_struc_data)
        # ------ output
    with open('cryspy.out', 'a') as fout:
        fout.write('{} structures by random\n'.format(rin.n_rand))

    # ---------- save init_struc_data
    pkl_data.save_init_struc(init_struc_data)

    # ----------return
    return init_struc_data, eagen
