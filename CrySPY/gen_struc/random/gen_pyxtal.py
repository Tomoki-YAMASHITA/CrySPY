'''
Random structure generation using PyXtal(https://github.com/qzhu2017/PyXtal)
'''

from contextlib import redirect_stdout
from multiprocessing import Process, Queue
import os
import random
import sys

import numpy as np
from pymatgen import DummySpecie, Structure, Molecule
from pyxtal import pyxtal
from pyxtal.database.collection import Collection

from ..struc_util import check_distance
from ..struc_util import sort_by_atype
from ..struc_util import out_poscar
from ..struc_util import scale_cell_mol
from ..struc_util import rot_mat


class Rnd_struc_gen_pyxtal:
    '''
    Random structure generation using pyxtal

    # ---------- args
    natot (int): number of atoms, e.g. 12 for Si4O8

    atype (list): atom type, e.g. ['Si', 'O'] for Si4O8

    nat (list): number of atom, e.g. [4, 8] for Si4O8

    vol_factor (list): default --> [1.0, 1.0]
                       [min, max] int or float

    vol_mu (int or float or None): default --> None
                                   average volume in Gaussian distribution
                                   when you scale cell volume

    vol_sigma (int or float or None): default --> None
                                      standard deviation in Gaussian distribution
                                      when you scale cell volume

    mindist (2d list): default --> None
                       constraint on minimum interatomic distance,
                       mindist must be a symmetric matrix
        e.g. [[1.8, 1.2], [1.2, 1.5]]
            Si - Si: 1.8 angstrom
            Si -  O: 1.2
             O -  O: 1.5

    spgnum ('all' or list): space group numbers which you use
                            'all' --> 1--230
                            list --> e.g. [1, 3, 100 ..., 229, 230]

    symprec (float): default --> 0.01
                     tolerance for symmetry finding
    '''

    def __init__(self, natot, atype, nat, vol_factor=[1.1, 1.1],
                 vol_mu=None, vol_sigma=None,
                 mindist=None, spgnum='all', symprec=0.01):
        # ---------- check args
        # ------ int
        for i in [natot]:
            if type(i) is int and i > 0:
                pass
            else:
                raise ValueError('natot must be positive int')
        # ------ list
        for x in [atype, nat, vol_factor]:
            if type(x) is not list:
                raise ValueError('atype, nat, and vol_factor must be list')
        if not len(atype) == len(nat):
            raise ValueError('not len(atype) == len(nat)')
        # ------ vol_factor
        if len(vol_factor) == 2:
            if not 0 < vol_factor[0] <= vol_factor[1]:
                raise ValueError('0 < vol_factor[0] <= vol_factor[1]')
        else:
            raise ValueError('length of vol_factor must be 2')
        # ------ vol_mu, vol_sigma
        if vol_mu is not None:
            for x in [vol_mu, vol_sigma]:
                if type(x) is not float and type(x) is not int:
                    raise ValueError('vol_mu and vol_sigma must be int or float')
            if vol_mu <= 0:
                raise ValueError('vol_mu must be positive')
        # ------ mindist
        if mindist is not None:
            if not len(atype) == len(mindist):
                raise ValueError('not len(atype) == len(mindist)')
            # -- check symmetric
            for i in range(len(mindist)):
                for j in range(i):
                    if not mindist[i][j] == mindist[j][i]:
                        raise ValueError('mindist is not symmetric. '
                                         '({}, {}): {}, ({}, {}): {}'.format(
                                             i, j, mindist[i][j],
                                             j, i, mindist[j][i]))

        # ------ spgnum
        if spgnum == 'all' or type(spgnum) is list:
            pass
        else:
            raise ValueError('spgnum is wrong. spgnum = {}'.format(spgnum))
        # ------ self.xxx = xxx
        self.natot = natot
        self.atype = atype
        self.nat = nat
        self.vol_factor = vol_factor
        self.vol_mu = vol_mu
        self.vol_sigma = vol_sigma
        self.mindist = mindist
        self.spgnum = spgnum
        self.symprec = symprec
        # ------ error list for spg number
        self.spg_error = []

    def set_mol(self, mol_file, nmol):
        '''
        set molecule files and number of molecules

        # ---------- args
        mol_file (list): list of path for molecular files (.xyz, etc)
                         one can also use pre-defined strings such as
                         H2O, CH4, benzene, ... etc.
                         see PyXtal document

        nmol (list): number of molecules

        # ---------- example
        # ------ 4 benzene molecules in unit cell
        mol_file = ['benzene']
        nmol = [4]
        # ------ molecules you make (2 * mol_1 and 2 * mol_2)
        mol_file = ['./mol_1.xyz', './mol_2.xyz']
        nmol = [2, 2]
        '''
        # ---------- check args
        # ------ mol_file
        if type(mol_file) is not list:
            raise ValueError('mol_file must be list')
        for s in mol_file:
            if type(s) is not str:
                raise ValueError('elements in mol_file must be string')
        # ------ nmol
        if type(nmol) is not list:
            raise ValueError('nmol must be list')
        for i in nmol:
            if type(i) is not int:
                raise ValueError('elements in nmol must be int')
        # ------ mol_file and nmol
        if not len(mol_file) == len(nmol):
            raise ValueError('not len(mol_file) == len(nmol)')
        # ----------
        mol_data = []
        pyxtal_mol_data = Collection('molecules')
        pyxtal_mol_names = list(Collection('molecules'))
        for i, mf in enumerate(mol_file):
            if os.path.isfile(mf):
                mol = Molecule.from_file(mf)
            elif mf in pyxtal_mol_names:
                mol = pyxtal_mol_data[mf]
            else:
                raise ValueError('no molecular files')
            mol_data.append(mol)
        # ---------- self.xxx
        self.nmol = nmol
        self.mol_data = mol_data

    def gen_struc(self, nstruc, id_offset=0, init_pos_path=None):
        '''
        Generate random structures for given space groups

        # ---------- args
        nstruc (int): number of generated structures

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        init_pos_path (str): default: None
                             specify a path of file
                             if you write POSCAR data of init_struc_data
                             ATTENSION: data are appended to the specified file

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- check args
        if not (type(nstruc) is int and nstruc > 0):
            raise ValueError('nstruc must be positive int')
        if type(id_offset) is not int:
            raise TypeError('id_offset must be int')
        if init_pos_path is None or type(init_pos_path) is str:
            pass
        else:
            raise ValueError('init_pos_path is wrong.'
                             ' init_pos_path = {}'.format(init_pos_path))
        # ---------- initialize
        self.init_struc_data = {}
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if self.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(self.spgnum)
            if spg in self.spg_error:
                continue
            # ------ vol_factor
            rand_vol = random.uniform(self.vol_factor[0], self.vol_factor[1])
            # ------ generate structure
            tmp_crystal = pyxtal()
            try:
                tmp_crystal.from_random(dim=3, group=spg, species=self.atype,
                                        numIons=self.nat, factor=rand_vol,
                                        conventional=False)
            except RuntimeError:
                sys.stderr.write('RuntimeError in spg = {} retry.\n'.format(spg))
                self.spg_error.append(spg)
                continue
            if tmp_crystal.valid:
                tmp_struc = tmp_crystal.to_pymatgen(resort=False)    # pymatgen Structure format
                # -- check the number of atoms
                if not self._check_nat(tmp_struc):
                    # (pyxtal 0.1.4) cryspy adopts "conventional=False",
                    #     which is better for DFT calculation
                    # pyxtal returns conventional cell, that is, too many atoms
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    if not self._check_nat(tmp_struc):    # failure
                        continue
                # -- sort, just in case
                tmp_struc = sort_by_atype(tmp_struc, self.atype)
                # -- scale volume
                if self.vol_mu is not None:
                    vol = random.gauss(mu=self.vol_mu, sigma=self.vol_sigma)
                    tmp_struc.scale_lattice(volume=vol)
                # -- check minimum distance
                if self.mindist is not None:
                    success, mindist_ij, dist = check_distance(tmp_struc,
                                                               self.atype,
                                                               self.mindist)
                    if not success:
                        sys.stderr.write('mindist in gen_struc: {} - {}, {}. retry.\n'.format(
                            self.atype[mindist_ij[0]],
                            self.atype[mindist_ij[1]],
                            dist))
                        continue    # failure
                # -- check actual space group
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=self.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                print('Structure ID {0:>6} was generated.'
                      ' Space group: {1:>3} --> {2:>3} {3}'.format(
                       cid, spg, spg_num, spg_sym))
                # -- save init_POSCARS
                if init_pos_path is not None:
                    out_poscar(tmp_struc, cid, init_pos_path)
            else:
                self.spg_error.append(spg)

    def gen_struc_mol(self, nstruc, id_offset=0, init_pos_path=None,
                      timeout_mol=180):
        '''
        Generate random molecular crystal structures for given space groups
        one have to run self.set_mol() in advance
        # ---------- args
        nstruc (int): number of generated structures

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        init_pos_path (str): default: None
                             specify a path of file
                             if you write POSCAR data of init_struc_data
                             ATTENSION: data are appended to the specified file

        timeout_mol (int or float): timeout for one molecular structure generation

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- check args
        if not (type(nstruc) is int and nstruc > 0):
            raise ValueError('nstruc must be positive int')
        if type(id_offset) is not int:
            raise TypeError('id_offset must be int')
        if init_pos_path is None or type(init_pos_path) is str:
            pass
        else:
            raise ValueError('init_pos_path is wrong.'
                             ' init_pos_path = {}'.format(init_pos_path))
        if type(timeout_mol) is not float and type(timeout_mol) is not int:
            raise ValueError('timeout_mol must be int or float')
        if timeout_mol <= 0:
            raise ValueError('timeout_mol must be positive')
        # ---------- initialize
        self.init_struc_data = {}
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if self.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(self.spgnum)
            if spg in self.spg_error:
                continue
            rand_vol = random.uniform(self.vol_factor[0], self.vol_factor[1])
            # ------ generate structure
            # -- multiprocess for measures against hangup
            q = Queue()
            p = Process(target=self._mp_mc, args=(spg, rand_vol, q))
            p.start()
            p.join(timeout=timeout_mol)
            if p.is_alive():
                p.terminate()
                p.join()
            if sys.version_info.minor >= 7:
                # Process.close() available from python 3.7
                p.close()
            if q.empty():
                sys.stderr.write('timeout for molecular structure generation. retry.\n')
                continue
            else:
                tmp_struc = q.get()
                tmp_valid = q.get()
                if tmp_struc == 'error':
                    self.spg_error.append(spg)
                    continue
            if tmp_valid:
                # -- scale volume
                if self.vol_mu is not None:
                    vol = random.gauss(mu=self.vol_mu, sigma=self.vol_sigma)
                    vol = vol * tmp_struc.num_sites / self.natot    # for conv. cell
                    tmp_struc = scale_cell_mol(tmp_struc, self.mol_data, vol)
                    if not tmp_struc:    # case: scale_cell_mol returns False
                        sys.stderr.write('failed scale cell. retry.\n')
                        continue
                # -- check nat
                if not self._check_nat(tmp_struc):
                    # cryspy adopts conventional=True
                    # pyxtal returns conventional cell,
                    # too many atoms if centering
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    if not self._check_nat(tmp_struc):    # failure
                        sys.stderr.write('different num. of atoms. retry.\n')
                        continue
                # -- sort, necessary in molecular crystal
                tmp_struc = sort_by_atype(tmp_struc, self.atype)
                # -- check minimum distance
                if self.mindist is not None:
                    success, mindist_ij, dist = check_distance(tmp_struc,
                                                               self.atype,
                                                               self.mindist)
                    if not success:
                        sys.stderr.write('mindist in gen_struc_mol: {} - {}, {}. retry.\n'.format(
                            self.atype[mindist_ij[0]],
                            self.atype[mindist_ij[1]],
                            dist))
                        continue    # failure
                # -- check actual space group
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=self.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                print('Structure ID {0:>6} was generated.'
                      ' Space group: {1:>3} --> {2:>3} {3}'.format(
                       cid, spg, spg_num, spg_sym))
                # -- save init_POSCARS
                if init_pos_path is not None:
                    out_poscar(tmp_struc, cid, init_pos_path)
            else:
                self.spg_error.append(spg)

    def gen_struc_mol_break_sym(self, nstruc, rot_mol=None,
                                id_offset=0, init_pos_path=None):
        '''
        Generate random molecular crystal structures
        one have to run self.set_mol() in advance
        molecules are put a Wyckoff position without consideration of symmetry

        # ---------- args
        nstruc (int): number of generated structures

        rot_mol (str): default: None
                       None, 'random', 'random_mol', or 'random_wyckoff'

        id_offset (int): default: 0
                         structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        init_pos_path (str): default: None
                             specify a path of file
                             if you write POSCAR data of init_struc_data
                             ATTENSION: data are appended to the specified file

        # ---------- comment
        generated structure data are saved in self.init_struc_data
        '''
        # ---------- check args
        if not (type(nstruc) is int and nstruc > 0):
            raise ValueError('nstruc must be positive int')
        if type(id_offset) is not int:
            raise TypeError('id_offset must be int')
        if init_pos_path is None or type(init_pos_path) is str:
            pass
        else:
            raise ValueError('init_pos_path is wrong.'
                             ' init_pos_path = {}'.format(init_pos_path))
        if rot_mol not in [None, 'random', 'random_mol', 'random_wyckoff']:
            raise ValueError('error in rot_mol')
        noble_gas = ['Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He']
        if len(self.mol_data) > len(noble_gas):
            raise ValueError('len(mol_data) > len(noble_gas)')
        # ---------- initialize
        self.init_struc_data = {}
        # ------ dummy atom type
        tmp_atype = noble_gas[:len(self.mol_data)]
        # ---------- loop for structure generattion
        while len(self.init_struc_data) < nstruc:
            # ------ spgnum --> spg
            if self.spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(self.spgnum)
            if spg in self.spg_error:
                continue
            # ------ vol_factor
            rand_vol = random.uniform(self.vol_factor[0], self.vol_factor[1])
            # ------ generate structure
            tmp_crystal = pyxtal()
            try:
                tmp_crystal.from_random(dim=3, group=spg, species=tmp_atype,
                                        numIons=self.nmol, factor=rand_vol,
                                        conventional=False)
            except RuntimeError:
                sys.stderr.write('RuntimeError in spg = {} retry.\n'.format(spg))
                self.spg_error.append(spg)
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
                if self.vol_mu is not None:
                    vol = random.gauss(mu=self.vol_mu, sigma=self.vol_sigma)
                    tmp_struc.scale_lattice(volume=vol)
                # -- rotate molecules
                if rot_mol == 'random_mol':
                    # each mol_data
                    mol_angles = []    # [angles of mol 1, angles of mol 2, ...]
                    for i in range(len(self.mol_data)):
                        mol_angles.append(2 * np.pi * np.random.rand(3))
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
                # -- replace dummy with mol
                dum_species = tmp_struc.species
                dum_coords = tmp_struc.cart_coords
                for (dum_specie, dum_coord) in zip(dum_species, dum_coords):
                    mol_index = tmp_atype.index(dum_type[dum_specie])
                    mol = self.mol_data[mol_index]
                    # rotation option
                    if rot_mol is None:
                        rot_mol_coord = mol.cart_coords
                    if rot_mol == 'random':
                        angle = 2 * np.pi * np.random.rand(3)
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    if rot_mol == 'random_mol':
                        angle = mol_angles[mol_index]
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    if rot_mol == 'random_wyckoff':
                        angle = dum_angles[dum_specie]
                        R = rot_mat(angle)
                        rot_mol_coord = np.matmul(mol.cart_coords, R)
                    # rotate coord
                    coord = rot_mol_coord + dum_coord
                    # append mol
                    for i, ms in enumerate(mol.species):
                        tmp_struc.append(ms, coord[i], coords_are_cartesian=True)
                # -- remove dummy
                tmp_struc.remove_sites(range(0, len(dum_species)))
                # -- check nat
                if not self._check_nat(tmp_struc):
                    # pyxtal returns conventional cell,
                    # too many atoms if centering
                    tmp_struc = tmp_struc.get_primitive_structure()
                    # recheck nat
                    if not self._check_nat(tmp_struc):    # failure
                        continue
                # -- sort, necessary in molecular crystal
                tmp_struc = sort_by_atype(tmp_struc, self.atype)
                # -- check minimum distance
                if self.mindist is not None:
                    success, mindist_ij, dist = check_distance(tmp_struc,
                                                               self.atype,
                                                               self.mindist)
                    if not success:
                        sys.stderr.write('mindist: {} - {}, {}. retry.\n'.format(
                            self.atype[mindist_ij[0]],
                            self.atype[mindist_ij[1]],
                            dist))
                        continue    # failure
                # -- check actual space group
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=self.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # -- register the structure in pymatgen format
                cid = len(self.init_struc_data) + id_offset
                self.init_struc_data[cid] = tmp_struc
                print('Structure ID {0:>6} was generated.'
                      ' Space group: {1:>3} --> {2:>3} {3}'.format(
                       cid, spg, spg_num, spg_sym))
                # -- save init_POSCARS
                if init_pos_path is not None:
                    out_poscar(tmp_struc, cid, init_pos_path)
            else:
                self.spg_error.append(spg)

    def _check_nat(self, struc):
        # ---------- count number of atoms in each element for check
        species_list = [a.species_string for a in struc]
        for i in range(len(self.atype)):
            if species_list.count(self.atype[i]) != self.nat[i]:
                return False    # failure
        return True

    def _mp_mc(self, spg, rand_vol, q):
        '''multiprocess part'''
        try:
            # ---------- temporarily stdout --> devnull
            with redirect_stdout(open(os.devnull, 'w')):
                tmp_crystal = pyxtal(molecular=True)
                tmp_crystal.from_random(dim=3, group=spg,
                                        species=self.mol_data, numIons=self.nmol,
                                        factor=rand_vol, conventional=False)
                # ---------- queue
            if tmp_crystal.valid:
                q.put(tmp_crystal.to_pymatgen(resort=False))
                q.put(tmp_crystal.valid)
            else:
                q.put(None)
                q.put(tmp_crystal.valid)
        except np.linalg.LinAlgError:
            sys.stderr.write('np.linalg.LinAlgError in spg = {} retry.\n'.format(spg))
            q.put('error')
            q.put(None)
        except RuntimeError:
            sys.stderr.write('RuntimeError in spg = {} retry.\n'.format(spg))
            q.put('error')
            q.put(None)
