'''
Random structure generation using PyXtal (https://github.com/qzhu2017/PyXtal)
'''
import collections
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from logging import getLogger
from multiprocessing import Process, Queue
import random
import sys

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie
from pyxtal import pyxtal
from pyxtal.tolerance import Tol_matrix

from ...util.struc_util import check_distance, sort_by_atype
from ...util.struc_util import get_nat, remove_zero, get_cn_comb
from ...util.struc_util import get_atype_dummy, scale_cell_mol, rot_mat


logger = getLogger('cryspy')


def gen_struc(
        nstruc,
        atype,
        nat,
        mindist,
        spgnum,
        symprec=0.01,
        id_offset=0,
        vol_factor=1.1,
        vol_mu=None,
        vol_sigma=None,
        vc=False,
        ll_nat=None,
        ul_nat=None,
        charge=None,
    ):
    '''
    Generate random structures for given space groups

        tuple may be replaced by list

    # ---------- args
    nstruc (int): number of structures to be generated
    atype (tuple): atom type (e.g. ('Na', 'Cl'))
    nat (tuple): number of atoms (e.g. (4, 4)), None if vc=True
    mindist (): minimum interatomic distance (e.g. ((2.0, 1.5), (1.5, 2.0)))
    spgnum (str, int, or tuple): space group number 'all', 0, or tuple of space group numbers
    symprec (float): symmetry tolerance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                             you obtain ID 10, ID 11, ID 12
    vol_factor (float): volume factor for structure generation
    vol_mu (float): mean for volume scaling
    vol_sigma (float): standard deviation for volume scaling
    # ------ for EA-vc
    vc (bool): variable composition. it needs ll_nat and ul_nat
    ll_nat (tuple): lower limit of number of atoms (e.g. (0, 0))
    ul_nat (tuple): upper limit of number of atoms (e.g. (8, 8))
    charge (tuple): charge of atoms (e.g. (1, -1)). Set if you want to check charge neutrality

    # ---------- return
    init_struc_data (dict): {ID: pymatgen Structure, ...}
    '''

    # ---------- initialize
    init_struc_data = {}
    if not vc:
        tmp_nat = nat
        tmp_atype = atype
        tolmat = _set_tol_mat(tmp_atype, mindist)
    if vc and charge is not None:
        cn_comb = get_cn_comb(ll_nat, ul_nat, charge)
        if not cn_comb:
            logger.error('No charge neutral combinations')
            raise SystemExit(1)
        logger.info(f'Consider charge neutrality: {charge}')

    # ---------- loop for structure generation
    while len(init_struc_data) < nstruc:
        # ------ spgnum --> spg
        if spgnum == 'all':
            spg = random.randint(1, 230)
        else:
            spg = random.choice(spgnum)
        # ------ generate structure
        tmp_crystal = pyxtal()
        if vc:    # variable composition
            if charge is None:
                nat = tuple([random.randint(l, u) for l, u in zip(ll_nat, ul_nat)])
            else:
                nat = random.choice(cn_comb)
            if sum(nat) == 0:
                continue    # restart
            if 0 in nat:    # remove 0 from numIons and corresponding index in atype, mindist
                tmp_atype, tmp_nat, tmp_mindist = remove_zero(atype, nat, mindist)
                tolmat = _set_tol_mat(tmp_atype, tmp_mindist)
            else:
                tmp_nat = tuple(nat)
                tmp_atype = atype
                tolmat = _set_tol_mat(tmp_atype, mindist)
        try:
            f = StringIO()    # to get output from pyxtal
            with redirect_stdout(f):
                tmp_crystal.from_random(
                    dim=3,
                    group=spg,
                    species=tmp_atype,
                    numIons=tmp_nat,
                    factor=vol_factor,
                    conventional=False,
                    tm=tolmat,
                )
            s = f.getvalue().rstrip()    # to delete \n
            if s:
                logger.warning(s)
        except Exception as e:
            logger.warning(f'{e}:     spg = {spg} retry.')
            continue
        if tmp_crystal.valid:
            tmp_struc = tmp_crystal.to_pymatgen(resort=False)
            tmp2_nat = get_nat(tmp_struc, tmp_atype)
            # -- check the number of atoms
            if tmp2_nat != tmp_nat:
                # (pyxtal 0.1.4) cryspy adopts "conventional=False",
                #     which is better for DFT calculation
                # pyxtal returns conventional cell, that is, too many atoms
                tmp_struc = tmp_struc.get_primitive_structure()
                # recheck nat
                tmp2_nat = get_nat(tmp_struc, tmp_atype)
                if tmp2_nat != tmp_nat:    # failure
                    continue
            # -- sort, just in case
            tmp_struc = sort_by_atype(tmp_struc, tmp_atype)
            # -- scale volume
            if vol_mu is not None:
                vol = random.gauss(mu=vol_mu, sigma=vol_sigma)
                tmp_struc.scale_lattice(volume=vol)
            # -- check actual space group
            try:
                spg_sym, spg_num = tmp_struc.get_space_group_info(
                    symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # -- tmp_struc --> init_struc_data
            cid = len(init_struc_data) + id_offset
            init_struc_data[cid] = tmp_struc
            logger.info(f'Structure ID {cid:>6}: {nat}'
                    f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    # ---------- return
    return init_struc_data


def gen_struc_mol(
        nstruc,
        atype,
        nat,
        mindist,
        spgnum,
        mol_data,
        nmol,
        symprec=0.01,
        id_offset=0,
        vol_factor=1.1,
        vol_mu=None,
        vol_sigma=None,
        timeout_mol=None,
    ):
    '''
    Generate random molecular crystal structures for given space groups

    # ---------- args
    nstruc (int): number of structures to be generated
    atype (tuple): atom type (e.g. ('Na', 'Cl'))
    nat (tuple): number of atoms (e.g. (4, 4))
    mindist (): minimum interatomic distance (e.g. ((2.0, 1.5), (1.5, 2.0)))
    spgnum (str, int, or tuple): space group number 'all', 0, or tuple of space group numbers
    mol_data (tuple or list): pymatgen Molecule data
    nmol (tuple): number of molecules (e.g. (4, 4))
    symprec (float): symmetry tolerance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                            you obtain ID 10, ID 11, ID 12
    vol_factor (float): volume factor for structure generation
    vol_mu (float): mean for volume scaling
    vol_sigma (float): standard deviation for volume scaling
    timeout_mol (None or float): if float, timeout for molecular structure generation

    # ---------- return
    init_struc_data (dict): {ID: pymatgen Structure, ...}
    struc_mol_id (dict):    not implemented yet, just return vacant dict for now
    '''

    # ---------- initialize
    init_struc_data = {}
    struc_mol_id = {}    # for EA and mol, not implemented yet

    # ---------- Tol_matrix
    tolmat = _set_tol_mat(atype, mindist)

    # ---------- loop for structure generattion
    while len(init_struc_data) < nstruc:
        # ------ spgnum --> spg
        if spgnum == 'all':
            spg = random.randint(1, 230)
        else:
            spg = random.choice(spgnum)
        # ------ generate structure
        if timeout_mol is None:
            tmp_crystal = pyxtal(molecular=True)
            try:
                f = StringIO()    # to get output from pyxtal
                with redirect_stdout(f):
                    tmp_crystal.from_random(
                        dim=3,
                        group=spg,
                        species=mol_data,
                        numIons=nmol,
                        factor=vol_factor,
                        conventional=False,
                        tm=tolmat
                    )
                s = f.getvalue().rstrip()    # to delete \n
                if s:
                    logger.warning(s)
                tmp_valid = tmp_crystal.valid
                if tmp_valid:
                    tmp_struc = tmp_crystal.to_pymatgen(resort=False)
            except Exception as e:
                logger.warning(f'{e}:     spg = {spg} retry.')
                continue
            # if algo in ['EA', 'EA-vc']:
            #     tmp_struc = tmp_crystal.to_pymatgen(resort=False)
            #     tmp_lattice = tmp_struc.lattice
            #     dums = []        # dummy atoms
            #     dum_pos = []     # internal position of dummy
            #     in_dists = []
            #     mol_dists = []
            #     for i, site in enumerate(tmp_crystal.mol_sites):
            #         for j, mol in enumerate(mol_data):
            #             if mol.species == site.mol.species:
            #                 dum = DummySpecie(f"X{j}{i}")
            #                 dums.append(dum)
            #                 dum_pos.append(site.position)
            #                 # -- calculate atom distance by pyxtal
            #                 in_dists.append(np.linalg.norm(site.position.dot(tmp_crystal.lattice.get_matrix()) -
            #                                                 site.get_coords_and_species(absolute=True)[0]
            #                                                 [0:len(site.mol.species)], axis=1))
            #     # -- calculate atom distance from self.mol_data
            #     for j, mol in enumerate(mol_data):
            #         tmp_dists = []
            #         for n, m in enumerate(mol):
            #             tmp_dists.append(mol.get_distance(0, n))
            #         mol_dists.append(tmp_dists)
            #     dum_struc = Structure.from_spacegroup(spg, tmp_lattice, dums, dum_pos)
            #     for ds, dc in zip(dum_struc.species, dum_struc.cart_coords):
            #         tmp_struc.append(ds, dc, coords_are_cartesian=True)
        else:
            # -- multiprocess for hangup prevention
            q = Queue()
            p = Process(target=_mp_mc, args=(tolmat, spg, mol_data, nmol,
                                            vol_factor, q))
            p.start()
            p.join(timeout=timeout_mol)
            if p.is_alive():
                p.terminate()
                p.join()
            if sys.version_info.minor >= 7:
                # Process.close() available from python 3.7
                p.close()
            if q.empty():
                logger.warning('timeout for molecular structure generation. retry.')
                continue
            else:
                # -- get struc data from _mp_mc
                # if rin.algo in ['EA', 'EA-vc']:
                #     tmp_q = q.get()
                #     tmp_struc = tmp_q[0]    # structure data
                #     dums = tmp_q[1]         # dummy element
                #     in_dists = tmp_q[2]     # interatomic distance in a molecule
                #     mol_dists = tmp_q[3]    # distance between molecules
                # else:
                #     tmp_struc = q.get()
                tmp_struc = q.get()
                tmp_valid = q.get()
                if tmp_struc == 'error':
                    # in case of 'error', tmp_valid <-- error message (Exception)
                    logger.warning(tmp_valid.args[0] + f': spg = {spg} retry.')
                    continue
        if tmp_valid:
            # -- scale volume
            if vol_mu is not None:
                vol = random.gauss(mu=vol_mu, sigma=vol_sigma)
                vol = vol * tmp_struc.num_sites / sum(nat)    # for conv. cell
                tmp_struc = scale_cell_mol(tmp_struc, mol_data, vol)
                if not tmp_struc:    # case: scale_cell_mol returns False
                    logger.warning('failed scale cell. retry.')
                    continue
            # -- check nat
            tmp_nat = get_nat(tmp_struc, atype)
            if tmp_nat != nat:
                # cryspy adopts conventional=True
                # pyxtal returns conventional cell,
                # too many atoms if centering
                tmp_struc = tmp_struc.get_primitive_structure()
                # recheck nat
                tmp_nat = get_nat(tmp_struc, atype)
                if tmp_nat != nat:    # failure
                    logger.warning(f'different num. of atoms. {tmp_nat}, {nat} retry.')
                    continue
            # -- grouping atoms for molecule using interatomic distance
            # if rin.algo in ['EA', 'EA-vc']:
            #     loopcnt = 0
            #     while loopcnt < 10:
            #         loopcnt += 1
            #         mol_group = []
            #         group_id = 0
            #         # check atom distance from dummy atom
            #         for i, ts in enumerate(tmp_struc.species):
            #             for dn, dummyatom in enumerate(dums):
            #                 if ts == dummyatom:
            #                     mol_indx = int(ts.symbol[1:2])
            #                     for j, tc_fc in enumerate(tmp_struc.frac_coords):
            #                         if i == j:
            #                             continue
            #                         dist = tmp_struc.get_distance(i, j)
            #                         for dist_indx, true_dist in enumerate(in_dists[dn]):
            #                             if round(true_dist, 2) == round(dist, 2):
            #                                 mol_group.append([tc_fc, group_id, mol_indx, dist_indx])
            #                                 break
            #                     group_id += 1
            #         success_cnt = 0
            #         # check number of atom per group
            #         # if failed, shuffle struc data, and continue(max 10times)
            #         mol_cnt = collections.Counter([x[2] for x in mol_group])
            #         for mc_key, mc_value in mol_cnt.items():
            #             if mc_value % len(mol_data[mc_key]) == 0:
            #                 success_cnt += 1
            #         if success_cnt == len(mol_data):
            #             break
            #         random.shuffle(tmp_struc)
            #     if loopcnt == 10:
            #         continue
            #     tmp_struc.remove_species(dums)
            # -- sort, necessary in molecular crystal
            tmp_struc = sort_by_atype(tmp_struc, atype)
            # -- record group_id, mol_id, and dist_index
            # if rin.algo in ['EA', 'EA-vc']:
            #     tmp_id = []
            #     tmp_mol_indx = []
            #     tmp_dist_indx = []
            #     for atom_fc in tmp_struc.frac_coords:
            #         for gatom in mol_group:
            #             if np.array_equal(atom_fc, gatom[0]):
            #                 tmp_id.append(gatom[1])
            #                 tmp_mol_indx.append(gatom[2])
            #                 tmp_dist_indx.append(gatom[3])
            # -- check minimum distance
            #    from CrySPY 0.10.4
            #    Tol_matrix is used for mindist
            #
            #if self.mindist is not None:
            #    success, mindist_ij, dist = check_distance(tmp_struc,
            #                                               self.atype,
            #                                               self.mindist)
            #    if not success:
            #        print('mindist in gen_struc_mol: {} - {}, {}. retry.'.format(
            #            self.atype[mindist_ij[0]],
            #            self.atype[mindist_ij[1]],
            #            dist), file=sys.stderr)
            #        continue    # failure
            # -- check actual space group
            try:
                spg_sym, spg_num = tmp_struc.get_space_group_info(
                    symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # -- register the structure in pymatgen format
            cid = len(init_struc_data) + id_offset
            init_struc_data[cid] = tmp_struc
            # if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
            #     struc_mol_id.update({cid: [tmp_mol_indx, tmp_id, mol_dists]})
            logger.info(f'Structure ID {cid:>6} was generated.'
                    f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    # ---------- return
    #            struc_mol_id is not implemented yet, just return vacant dict for now
    return init_struc_data, struc_mol_id


def gen_struc_mol_break_sym(
        nstruc,
        atype,
        nat,
        mindist,
        mindist_dummy,
        spgnum,
        mol_data,
        nmol,
        symprec=0.01,
        id_offset=0,
        vol_factor=1.1,
        vol_mu=None,
        vol_sigma=None,
        rot_mol='random_wyckoff',
        nrot=20,
    ):
    '''
    Generate random molecular crystal structures
    molecules are put at a Wyckoff position without consideration of symmetry

    # ---------- args
    nstruc (int): number of structures to be generated
    atype (tuple): atom type (e.g. ('Na', 'Cl'))
    nat (tuple): number of atoms (e.g. (4, 4))
    mindist (): minimum interatomic distance (e.g. ((2.0, 1.5), (1.5, 2.0)))
    mindist_dummy: mindist for dummy atoms (= intermolecular distance)
    spgnum (str, int, or tuple): space group number 'all', 0, or tuple of space group numbers
    mol_data (tuple or list): pymatgen Molecule data
    nmol (tuple): number of molecules (e.g. (4, 4))
    symprec (float): symmetry tolerance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                            you obtain ID 10, ID 11, ID 12
    vol_factor (float): volume factor for structure generation
    vol_mu (float): mean for volume scaling
    vol_sigma (float): standard deviation for volume scaling
    rot_mol (str): rotation option. 'random', 'random_mol', or 'random_wyckoff'
    nrot (int): maximum number of trials to rotate molecules

    # ---------- return
    init_struc_data (dict): {ID: pymatgen Structure, ...}
    struc_mol_id (dict):   not implemented yet, just return vacant dict for now
    '''

    # ---------- initialize
    init_struc_data = {}
    struc_mol_id = {}    # for EA and mol, not implemented yet

    # ------ dummy atom type
    atype_dummy = get_atype_dummy(len(mol_data))

    # ---------- Tol_matrix for dummy atoms
    tolmat = _set_tol_mat(atype_dummy, mindist_dummy)

    # ---------- loop for structure generattion
    while len(init_struc_data) < nstruc:
        # ------ spgnum --> spg
        if spgnum == 'all':
            spg = random.randint(1, 230)
        else:
            spg = random.choice(spgnum)
        # ------ generate structure
        tmp_crystal = pyxtal()
        try:
            f = StringIO()
            with redirect_stdout(f):
                tmp_crystal.from_random(
                    dim=3,
                    group=spg,
                    species=atype_dummy,
                    numIons=nmol,
                    factor=vol_factor,
                    conventional=False,
                    tm=tolmat
                )
            s = f.getvalue().rstrip()    # to delete \n
            if s:
                logger.warning(s)
        except Exception as e:
            logger.warning(f'{e}:    spg = {spg} retry.')
            continue
        if tmp_crystal.valid:
            # -- each wyckoff position --> dummy atom
            dums = []        # dummy atoms
            dum_pos = []     # internal position of dummy
            dum_type = {}    # type of dummy
                                #  e.g. {DummySpecie X00+: 'Rn',
                                #        DummySpecie X10+: 'Rn',
                                #        DummySpecie X20+: 'Xe',
                                #        DummySpecie X30+: 'Xe'}
            for i, site in enumerate(tmp_crystal.atom_sites):
                dum = DummySpecie(f"X{i}")
                dums.append(dum)
                dum_pos.append(site.position)
                dum_type[dum] = site.specie
            # -- tmp_struc with dummy atoms
            lattice = tmp_crystal.lattice.get_matrix()
            tmp_struc = Structure.from_spacegroup(spg, lattice, dums, dum_pos)
            tmp_struc = tmp_struc.get_primitive_structure()
            # -- scale volume
            if vol_mu is not None:
                vol = random.gauss(mu=vol_mu, sigma=vol_sigma)
                tmp_struc.scale_lattice(volume=vol)
            # -- save dummy coords
            dum_species = tmp_struc.species
            dum_coords = tmp_struc.cart_coords
            # if rin.algo in ['EA', 'EA-vc']:
            #     tmp_mol_indx = []
            #     tmp_id = []
            #     tmp_id_cnt = 0
            # -- remove dummy
            tmp_struc.remove_sites(range(0, len(dum_species)))
            tmp_struc_ori = tmp_struc.copy()
            rot_success = False
            # -- rotate molecules
            for nrel in range(nrot):
                tmp_struc = tmp_struc_ori.copy()
                for (dum_specie, dum_coord) in zip(dum_species, dum_coords):
                    mol_index = atype_dummy.index(dum_type[dum_specie])
                    mol = mol_data[mol_index]
                    # rotation option
                    if rot_mol is None:
                        rot_mol_coord = mol.cart_coords
                    if rot_mol == 'random':
                        angle = 2 * np.pi * np.random.rand(3)
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    if rot_mol == 'random_mol':
                        # each mol_data
                        mol_angles = []    # [angles of mol 1, angles of mol 2, ...]
                        for i in range(len(mol_data)):
                            mol_angles.append(2 * np.pi * np.random.rand(3))
                        angle = mol_angles[mol_index]
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    if rot_mol == 'random_wyckoff':
                        # each Wyckoff
                        dum_angles = {}    # e.g.
                                            # {DummySpecie X00+: array([ , , ]),
                                            #  DummySpecie X10+: array([ , , ]),
                                            #  DummySpecie X20+: array([ , , ]),
                                            #  DummySpecie X30+: array([ , , ])}
                        angles = 2 * np.pi * np.random.rand(len(dums), 3)
                        for (dum, angle) in zip(dums, angles):
                            dum_angles[dum] = angle
                        angle = dum_angles[dum_specie]
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    # rotate coord
                    coord = rot_mol_coord + dum_coord
                    # append mol
                    for i, ms in enumerate(mol.species):
                        tmp_struc.append(ms, coord[i], coords_are_cartesian=True)
                        # if rin.algo in ['EA', 'EA-vc']:
                        #     tmp_id.append(tmp_id_cnt)
                        #     tmp_mol_indx.append(mol_index)
                    # if rin.algo in ['EA', 'EA-vc']:
                    #     tmp_id_cnt += 1
                # -- check nat
                tmp_nat = get_nat(tmp_struc, atype)
                if tmp_nat != nat:
                    # pyxtal returns conventional cell,
                    # too many atoms if centering
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    tmp_nat = get_nat(tmp_struc, atype)
                    if tmp_nat != nat:    # failure
                        if rot_mol is None:
                            break    # go back to the while loop
                        continue
                # -- sort, necessary in molecular crystal
                tmp_struc = sort_by_atype(tmp_struc, atype)
                # -- check minimum distance
                if mindist is not None:
                    success, mindist_ij, dist = check_distance(tmp_struc,
                                                               atype,
                                                               mindist)
                    if not success:
                        type0 = atype[mindist_ij[0]]
                        type1 = atype[mindist_ij[1]]
                        logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
                        if rot_mol is None:
                            break    # go back to the while loop
                        continue    # failure
                # -- check actual space group (success)
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                rot_success = True
                break
            # -- reach maximum times to rotate (failure)
            if not rot_success:
                continue    # go back to the while loop
            # -- tmp_struc --> init_struc_data
            cid = len(init_struc_data) + id_offset
            init_struc_data[cid] = tmp_struc
            # if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
            #     # calculate interatomic distance of mol_data
            #     mol_dists = []
            #     for j, mol in enumerate(mol_data):
            #         tmp_dists = []
            #         for n, m in enumerate(mol):
            #             tmp_dists.append(mol.get_distance(0, n))
            #         mol_dists.append(tmp_dists)
            #     struc_mol_id.update({cid: [tmp_mol_indx, tmp_id, mol_dists]})
            logger.info(f'Structure ID {cid:>6} was generated.'
                    f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    # ---------- return
    #            struc_mol_id is not implemented yet, just return vacant dict for now
    return init_struc_data, struc_mol_id


def _set_tol_mat(atype, mindist):
    tolmat = Tol_matrix()
    for i, itype in enumerate(atype):
        for j, jtype in enumerate(atype):
            if i <= j:
                tolmat.set_tol(itype, jtype, mindist[i][j])
    return tolmat


def _mp_mc(tolmat, spg, mol_data, nmol, vol_factor, q):
    '''
    multiprocess part
    here cannot use rin.xxx and logging
    '''
    try:
        np.random.seed(random.randint(0, 100000000))
        tmp_crystal = pyxtal(molecular=True)
        f = StringIO()
        with redirect_stdout(f):
            with redirect_stderr(f):
                tmp_crystal.from_random(dim=3, group=spg,
                                        species=mol_data, numIons=nmol,
                                        factor=vol_factor, conventional=False, tm=tolmat)
        s = f.getvalue().rstrip()    # to delete \n
        if s:
            logger.warning(s)
        # if algo in ['EA', 'EA-vc']:
        #     tmp_struc = tmp_crystal.to_pymatgen(resort=False)
        #     tmp_lattice = tmp_struc.lattice
        #     dums = []        # dummy atoms
        #     dum_pos = []     # internal position of dummy
        #     in_dists = []
        #     mol_dists = []
        #     for i, site in enumerate(tmp_crystal.mol_sites):
        #         for j, mol in enumerate(mol_data):
        #             if mol.species == site.mol.species:
        #                 dum = DummySpecie(f"X{j}{i}")
        #                 dums.append(dum)
        #                 dum_pos.append(site.position)
        #                 # -- calculate atom distance by pyxtal
        #                 in_dists.append(np.linalg.norm(site.position.dot(tmp_crystal.lattice.get_matrix()) -
        #                                                 site.get_coords_and_species(absolute=True)[0]
        #                                                 [0:len(site.mol.species)], axis=1))
        #     # -- calculate atom distance from self.mol_data
        #     for j, mol in enumerate(mol_data):
        #         tmp_dists = []
        #         for n, m in enumerate(mol):
        #             tmp_dists.append(mol.get_distance(0, n))
        #         mol_dists.append(tmp_dists)
        #     dum_struc = Structure.from_spacegroup(spg, tmp_lattice, dums, dum_pos)
        #     for ds, dc in zip(dum_struc.species, dum_struc.cart_coords):
        #         tmp_struc.append(ds, dc, coords_are_cartesian=True)
        # ---------- queue
        if tmp_crystal.valid:
            # if algo in ['EA', 'EA-vc']:
            #     q.put([tmp_struc, dums, in_dists, mol_dists])
            # else:
            #     q.put(tmp_crystal.to_pymatgen(resort=False))
            q.put(tmp_crystal.to_pymatgen(resort=False))
            q.put(tmp_crystal.valid)
        else:
            q.put(None)
            q.put(tmp_crystal.valid)
    except Exception as e:
        q.put('error')
        q.put(e)
