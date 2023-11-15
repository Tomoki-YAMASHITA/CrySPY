'''
Random structure generation using PyXtal (https://github.com/qzhu2017/PyXtal)
'''

import collections
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
from logging import getLogger
from multiprocessing import Process, Queue
import os
import random
import sys

import numpy as np
from pymatgen.core import Structure, Molecule
from pymatgen.core.periodic_table import DummySpecie
from pyxtal import pyxtal
from pyxtal.database.collection import Collection
from pyxtal.tolerance import Tol_matrix

from ...util.struc_util import check_distance, sort_by_atype, get_nat
from ...util.struc_util import get_atype_dummy, scale_cell_mol, rot_mat
from ...IO import read_input as rin


logger = getLogger('cryspy')

class Rnd_struc_gen_pyxtal:
    '''
    Random structure generation using pyxtal
    '''

    def __init__(self, mindist):
        self.mindist = mindist

    def set_mol(self):
        '''
        set molecule files and number of molecules
        '''
        # ----------
        mol_data = []
        pyxtal_mol_data = Collection('molecules')
        pyxtal_mol_names = list(Collection('molecules'))
        for i, mf in enumerate(rin.mol_file):
            if os.path.isfile(mf):
                mol = Molecule.from_file(mf)
            elif mf in pyxtal_mol_names:
                mol = pyxtal_mol_data[mf]
            else:
                logger.error('no molecular files')
                raise SystemExit(1)
            mol_data.append(mol)
        # ---------- self.xxx
        self.mol_data = mol_data

    def gen_struc(self, nstruc, id_offset=0, vc=False):
        '''
        Generate random structures for given space groups

        # ---------- args
        nstruc (int): number of generated structures

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        vc (bool): variable composition. it needs ll_nat and ul_nat

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- initialize
        self.init_struc_data = {}
        # ---------- Tol_matrix
        tolmat = self._set_tol_mat(rin.atype, self.mindist)
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if rin.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(rin.spgnum)
            # ------ vol_factor
            rand_vol = random.uniform(rin.vol_factor[0], rin.vol_factor[1])
            # ------ generate structure
            tmp_crystal = pyxtal()
            if not vc:
                numIons = rin.nat
            else:    # variable composition
                numIons = [random.randint(l, u) for l, u in zip(rin.ll_nat, rin.ul_nat)]
            try:
                f = StringIO()
                with redirect_stdout(f):
                    tmp_crystal.from_random(dim=3, group=spg, species=rin.atype,
                                            numIons=numIons, factor=rand_vol,
                                            conventional=False, tm=tolmat)
                s = f.getvalue().rstrip()    # to delete \n
                if s:
                    logger.warning(s)
            except Exception as e:
                logger.warning(str(e.args[0]) + f': spg = {spg} retry.')
                continue
            if tmp_crystal.valid:
                tmp_struc = tmp_crystal.to_pymatgen(resort=False)    # pymatgen Structure format
                tmp_nat, _ = get_nat(tmp_struc, rin.atype)
                # -- check the number of atoms
                if tmp_nat != numIons:
                    # (pyxtal 0.1.4) cryspy adopts "conventional=False",
                    #     which is better for DFT calculation
                    # pyxtal returns conventional cell, that is, too many atoms
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    if tmp_nat != numIons:    # failure
                        continue
                # -- sort, just in case
                tmp_struc = sort_by_atype(tmp_struc, rin.atype)
                # -- scale volume
                if rin.vol_mu is not None:
                    vol = random.gauss(mu=rin.vol_mu, sigma=rin.vol_sigma)
                    tmp_struc.scale_lattice(volume=vol)
                # -- check minimum distance
                #    from CrySPY 0.10.4
                #    Tol_matrix is used for mindist
                #
                #if self.mindist is not None:
                #    success, mindist_ij, dist = check_distance(tmp_struc,
                #                                               self.atype,
                #                                               self.mindist)
                #    if not success:
                #        print('mindist in gen_struc: {} - {}, {}. retry.'.format(
                #            self.atype[mindist_ij[0]],
                #            self.atype[mindist_ij[1]],
                #            dist), file=sys.stderr)
                #        continue    # failure
                # -- check actual space group
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                logger.info(f'Structure ID {cid:>6}: {numIons}'
                      f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    def gen_struc_mol(self, nstruc, id_offset=0):
        '''
        Generate random molecular crystal structures for given space groups
        one have to run self.set_mol() in advance
        # ---------- args
        nstruc (int): number of generated structures

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- initialize
        self.init_struc_data = {}
        if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
            self.struc_mol_id = {}
        # ---------- Tol_matrix
        tolmat = self._set_tol_mat(rin.atype, self.mindist)
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if rin.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(rin.spgnum)
            rand_vol = random.uniform(rin.vol_factor[0], rin.vol_factor[1])
            # ------ generate structure
            # -- multiprocess for measures against hangup
            q = Queue()
            p = Process(target=self._mp_mc, args=(tolmat, spg, rin.nmol, rand_vol, q, rin.algo))
            p.start()
            p.join(timeout=rin.timeout_mol)
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
                if rin.algo in ['EA', 'EA-vc']:
                    tmp_q = q.get()
                    tmp_struc = tmp_q[0]    # structure data
                    dums = tmp_q[1]         # dummy element
                    in_dists = tmp_q[2]     # interatomic distance in a molecule
                    mol_dists = tmp_q[3]    # distance between molecules
                else:
                    tmp_struc = q.get()
                tmp_valid = q.get()
                if tmp_struc == 'error':
                    # in case of 'error', tmp_valid <-- error message (Exception)
                    logger.warning(tmp_valid.args[0] + f': spg = {spg} retry.')
                    continue
            if tmp_valid:
                # -- scale volume
                if rin.vol_mu is not None:
                    vol = random.gauss(mu=rin.vol_mu, sigma=rin.vol_sigma)
                    vol = vol * tmp_struc.num_sites / rin.natot    # for conv. cell
                    tmp_struc = scale_cell_mol(tmp_struc, self.mol_data, vol)
                    if not tmp_struc:    # case: scale_cell_mol returns False
                        logger.warning('failed scale cell. retry.')
                        continue
                # -- check nat
                tmp_nat, _ = get_nat(tmp_struc, rin.atype)
                if tmp_nat != rin.nat:
                    # cryspy adopts conventional=True
                    # pyxtal returns conventional cell,
                    # too many atoms if centering
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    if tmp_nat != rin.nat:    # failure
                        logger.warning('different num. of atoms. retry.')
                        continue
                # -- grouping atoms for molecule using interatomic distance
                if rin.algo in ['EA', 'EA-vc']:
                    loopcnt = 0
                    while loopcnt < 10:
                        loopcnt += 1
                        mol_group = []
                        group_id = 0
                        # check atom distance from dummy atom
                        for i, ts in enumerate(tmp_struc.species):
                            for dn, dummyatom in enumerate(dums):
                                if ts == dummyatom:
                                    mol_indx = int(ts.symbol[1:2])
                                    for j, tc_fc in enumerate(tmp_struc.frac_coords):
                                        if i == j:
                                            continue
                                        dist = tmp_struc.get_distance(i, j)
                                        for dist_indx, true_dist in enumerate(in_dists[dn]):
                                            if round(true_dist, 2) == round(dist, 2):
                                                mol_group.append([tc_fc, group_id, mol_indx, dist_indx])
                                                break
                                    group_id += 1
                        success_cnt = 0
                        # check number of atom per group
                        # if failed, shuffle struc data, and continue(max 10times)
                        mol_cnt = collections.Counter([x[2] for x in mol_group])
                        for mc_key, mc_value in mol_cnt.items():
                            if mc_value % len(self.mol_data[mc_key]) == 0:
                                success_cnt += 1
                        if success_cnt == len(self.mol_data):
                            break
                        random.shuffle(tmp_struc)
                    if loopcnt == 10:
                        continue
                    tmp_struc.remove_species(dums)
                # -- sort, necessary in molecular crystal
                tmp_struc = sort_by_atype(tmp_struc, rin.atype)
                # -- record group_id, mol_id, and dist_index
                if rin.algo in ['EA', 'EA-vc']:
                    tmp_id = []
                    tmp_mol_indx = []
                    tmp_dist_indx = []
                    for atom_fc in tmp_struc.frac_coords:
                        for gatom in mol_group:
                            if np.array_equal(atom_fc, gatom[0]):
                                tmp_id.append(gatom[1])
                                tmp_mol_indx.append(gatom[2])
                                tmp_dist_indx.append(gatom[3])
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
                        symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
                    self.struc_mol_id.update({cid: [tmp_mol_indx, tmp_id, mol_dists]})
                logger.info(f'Structure ID {cid:>6} was generated.'
                      f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    def gen_struc_mol_break_sym(self, nstruc, mindist_dummy,
                                id_offset=0):
        '''
        Generate random molecular crystal structures
        one have to run self.set_mol() in advance
        molecules are put a Wyckoff position without consideration of symmetry

        # ---------- args
        nstruc (int): number of generated structures

        mindist_dummy: mindist for dummy atoms

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- initialize
        self.init_struc_data = {}
        if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
            self.struc_mol_id = {}
        # ------ dummy atom type
        atype_dummy = get_atype_dummy()
        # ---------- Tol_matrix for dummy atoms
        tolmat = self._set_tol_mat(atype_dummy, mindist_dummy)
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if rin.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(rin.spgnum)
            # ------ vol_factor
            rand_vol = random.uniform(rin.vol_factor[0], rin.vol_factor[1])
            # ------ generate structure
            tmp_crystal = pyxtal()
            try:
                f = StringIO()
                with redirect_stdout(f):
                    tmp_crystal.from_random(dim=3, group=spg, species=atype_dummy,
                                            numIons=rin.nmol, factor=rand_vol,
                                            conventional=False, tm=tolmat)
                s = f.getvalue().rstrip()    # to delete \n
                if s:
                    logger.warning(s)
            except Exception as e:
                logger.warning(str(e.args[0]) + f': spg = {spg} retry.')
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
                    dum = DummySpecie("X{}".format(i))
                    dums.append(dum)
                    dum_pos.append(site.position)
                    dum_type[dum] = site.specie
                # -- tmp_struc with dummy atoms
                lattice = tmp_crystal.lattice.get_matrix()
                tmp_struc = Structure.from_spacegroup(spg, lattice, dums, dum_pos)
                tmp_struc = tmp_struc.get_primitive_structure()
                # -- scale volume
                if rin.vol_mu is not None:
                    vol = random.gauss(mu=rin.vol_mu, sigma=rin.vol_sigma)
                    tmp_struc.scale_lattice(volume=vol)
                # -- save dummy coords
                dum_species = tmp_struc.species
                dum_coords = tmp_struc.cart_coords
                if rin.algo in ['EA', 'EA-vc']:
                    tmp_mol_indx = []
                    tmp_id = []
                    tmp_id_cnt = 0
                # -- remove dummy
                tmp_struc.remove_sites(range(0, len(dum_species)))
                tmp_struc_ori = tmp_struc.copy()
                rot_success = False
                # -- rotate molecules
                for nrel in range(rin.nrot):
                    tmp_struc = tmp_struc_ori.copy()
                    for (dum_specie, dum_coord) in zip(dum_species, dum_coords):
                        mol_index = atype_dummy.index(dum_type[dum_specie])
                        mol = self.mol_data[mol_index]
                        # rotation option
                        if rin.rot_mol is None:
                            rot_mol_coord = mol.cart_coords
                        if rin.rot_mol == 'random':
                            angle = 2 * np.pi * np.random.rand(3)
                            R = rot_mat(angle)
                            rot_mol_coord = np.matmul(mol.cart_coords, R)
                        if rin.rot_mol == 'random_mol':
                            # each mol_data
                            mol_angles = []    # [angles of mol 1, angles of mol 2, ...]
                            for i in range(len(self.mol_data)):
                                mol_angles.append(2 * np.pi * np.random.rand(3))
                            angle = mol_angles[mol_index]
                            R = rot_mat(angle)
                            rot_mol_coord = np.matmul(mol.cart_coords, R)
                        if rin.rot_mol == 'random_wyckoff':
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
                            if rin.algo in ['EA', 'EA-vc']:
                                tmp_id.append(tmp_id_cnt)
                                tmp_mol_indx.append(mol_index)
                        if rin.algo in ['EA', 'EA-vc']:
                            tmp_id_cnt += 1
                    # -- check nat
                    tmp_nat, _ = get_nat(tmp_struc, rin.atype)
                    if tmp_nat != rin.nat:
                        # pyxtal returns conventional cell,
                        # too many atoms if centering
                        tmp_struc = tmp_struc.get_primitive_structure()
                        # recheck nat
                        if tmp_nat != rin.nat:    # failure
                            if rin.rot_mol is None:
                                break    # go back to the while loop
                            continue
                    # -- sort, necessary in molecular crystal
                    tmp_struc = sort_by_atype(tmp_struc, rin.atype)
                    # -- check minimum distance
                    if self.mindist is not None:
                        success, mindist_ij, dist = check_distance(tmp_struc,
                                                                   rin.atype,
                                                                   self.mindist)
                        if not success:
                            type0 = rin.atype[mindist_ij[0]]
                            type1 = rin.atype[mindist_ij[1]]
                            logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
                            if rin.rot_mol is None:
                                break    # go back to the while loop
                            continue    # failure
                    # -- check actual space group (success)
                    try:
                        spg_sym, spg_num = tmp_struc.get_space_group_info(
                            symprec=rin.symprec)
                    except TypeError:
                        spg_num = 0
                        spg_sym = None
                    rot_success = True
                    break
                # -- reach maximum times to rotate (failure)
                if not rot_success:
                    continue    # go back to the while loop
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
                    # calculate interatomic distance of mol_data
                    mol_dists = []
                    for j, mol in enumerate(self.mol_data):
                        tmp_dists = []
                        for n, m in enumerate(mol):
                            tmp_dists.append(mol.get_distance(0, n))
                        mol_dists.append(tmp_dists)
                    self.struc_mol_id.update({cid: [tmp_mol_indx, tmp_id, mol_dists]})
                logger.info(f'Structure ID {cid:>6} was generated.'
                      f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')

    def _set_tol_mat(self, atype, mindist):
        tolmat = Tol_matrix()
        for i, itype in enumerate(atype):
            for j, jtype in enumerate(atype):
                if i <= j:
                    tolmat.set_tol(itype, jtype, mindist[i][j])
        return tolmat

    def _mp_mc(self, tolmat, spg, nmol, rand_vol, q, algo):
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
                                            species=self.mol_data, numIons=nmol,
                                            factor=rand_vol, conventional=False, tm=tolmat)
            s = f.getvalue().rstrip()    # to delete \n
            if s:
                logger.warning(s)
            if algo in ['EA', 'EA-vc']:
                tmp_struc = tmp_crystal.to_pymatgen(resort=False)
                tmp_lattice = tmp_struc.lattice
                dums = []        # dummy atoms
                dum_pos = []     # internal position of dummy
                in_dists = []
                mol_dists = []
                for i, site in enumerate(tmp_crystal.mol_sites):
                    for j, mol in enumerate(self.mol_data):
                        if mol.species == site.mol.species:
                            dum = DummySpecie("X{}{}".format(j, i))
                            dums.append(dum)
                            dum_pos.append(site.position)
                            # -- calculate atom distance by pyxtal
                            in_dists.append(np.linalg.norm(site.position.dot(tmp_crystal.lattice.get_matrix()) -
                                                           site.get_coords_and_species(absolute=True)[0]
                                                           [0:len(site.mol.species)], axis=1))
                # -- calculate atom distance from self.mol_data
                for j, mol in enumerate(self.mol_data):
                    tmp_dists = []
                    for n, m in enumerate(mol):
                        tmp_dists.append(mol.get_distance(0, n))
                    mol_dists.append(tmp_dists)
                dum_struc = Structure.from_spacegroup(spg, tmp_lattice, dums, dum_pos)
                for ds, dc in zip(dum_struc.species, dum_struc.cart_coords):
                    tmp_struc.append(ds, dc, coords_are_cartesian=True)
            # ---------- queue
            if tmp_crystal.valid:
                if algo in ['EA', 'EA-vc']:
                    q.put([tmp_struc, dums, in_dists, mol_dists])
                else:
                    q.put(tmp_crystal.to_pymatgen(resort=False))
                q.put(tmp_crystal.valid)
            else:
                q.put(None)
                q.put(tmp_crystal.valid)
        except Exception as e:
            q.put('error')
            q.put(e)
