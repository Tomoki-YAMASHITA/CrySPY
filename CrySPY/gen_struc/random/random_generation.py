'''
Random structure generation
'''

import json
import math
import os
import random
import subprocess
import sys

import numpy as np
from pymatgen import Structure

from ..struc_util import check_distance
from ..struc_util import out_poscar


class Rnd_struc_gen:
    '''
    Random structure generation

    # ---------- args
    natot (int): number of atoms, e.g. 12 for Si4O8

    atype (list): atom type, e.g. ['Si', 'O'] for Si4O8

    nat (list): number of atom, e.g. [4, 8] for Si4O8

    minlen (float): minimum length of lattice vector

    maxlen (float): maximum length of lattice vector

    dangle (float): delta angle for alpha, beta, and gamma in degree unit

    mindist (2d list): constraint on minimum interatomic distance,
                       mindist must be a symmetric matrix
        e.g. [[1.8, 1.2], [1.2, 1.5]]
            Si - Si: 1.8 angstrom
            Si -  O: 1.2
             O -  O: 1.5

    vol_mu (int or float or None): default --> None
                                   average volume in Gaussian distribution
                                   when you scale cell volume

    vol_sigma (int or float or None): default --> None
                                      standard deviation in Gaussian distribution
                                      when you scale cell volume

    symprec (float): default --> 0.01
                     tolerance for symmetry finding

    # ---------- instance methods
    self.gen_wo_spg(self, nstruc)

    self.gen_with_find_wy(self, nstruc, spgnum='all')
    '''

    def __init__(self, natot, atype, nat, minlen, maxlen, dangle, mindist,
                 vol_mu=None, vol_sigma=None, maxcnt=50, symprec=0.01):
        # ---------- check args
        # ------ int
        for i in [natot, maxcnt]:
            if type(i) is int and i > 0:
                pass
            else:
                raise ValueError('natot and maxcnt must be positive int')
        # ------ list
        for x in [atype, nat, mindist]:
            if type(x) is not list:
                raise ValueError('atype, nat, and mindist must be list')
        if not len(atype) == len(nat) == len(mindist):
            raise ValueError('not len(atype) == len(nat) == len(mindist)')
        # -- check symmetric
        for i in range(len(mindist)):
            for j in range(i):
                if not mindist[i][j] == mindist[j][i]:
                    raise ValueError('mindist is not symmetric. '
                                     '({}, {}): {}, ({}, {}): {}'.format(
                                         i, j, mindist[i][j],
                                         j, i, mindist[j][i]))
        # ------ vol_mu, vol_sigma
        if vol_mu is not None:
            for x in [vol_mu, vol_sigma]:
                if type(x) is not float and type(x) is not int:
                    raise ValueError('vol_mu and vol_sigma must be int or float')
            if vol_mu <= 0:
                raise ValueError('vol_mu must be positive')
        # ------ float
        for x in [minlen, maxlen, dangle, symprec]:
            if minlen > maxlen:
                raise ValueError('minlen > maxlen')
            if type(x) is int:
                x = float(x)    # int --> float
            if type(x) is float and x > 0.0:
                pass
            else:
                raise ValueError('minlen, maxlen, dangle, and symprec'
                                 ' must be positive float')
        # ------ self.xxx = xxx
        self.natot = natot
        self.maxcnt = maxcnt
        self.atype = atype
        self.nat = nat
        self.mindist = mindist
        self.vol_mu = vol_mu
        self.vol_sigma = vol_sigma
        self.minlen = minlen
        self.maxlen = maxlen
        self.dangle = dangle
        self.symprec = symprec

    def gen_wo_spg(self, nstruc, id_offset=0, init_pos_path=None):
        '''
        Generate random structures without space group information

        # ---------- args
        nstruc (int): number of generated structures

        id_offset (int): structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        init_pos_path (str): specify a path of file,
                             if you write POSCAR data of init_struc_data
                             ATTENSION: data are appended to the specified file

        # ---------- comment
        generated init_struc_data is saved in self.init_struc_data
        '''
        # ---------- check args
        if type(nstruc) is int and nstruc > 0:
            pass
        else:
            raise ValueError('nstruc must be positive int')
        if type(id_offset) is not int:
            raise TypeError('id_offset must be int')
        if init_pos_path is None or type(init_pos_path) is str:
            pass
        else:
            raise ValueError('init_pos_path is wrong.'
                             ' init_pos_path = {}'.format(init_pos_path))
        # ---------- initialize
        init_struc_data = {}
        self._get_atomlist()    # get self.atomlist
        # ---------- generate structures
        while len(init_struc_data) < nstruc:
            # ------ get spg, a, b, c, alpha, beta, gamma in self.*
            self._gen_lattice(spgnum=0)
            # ------ get va, vb, and vc in self.*
            self._calc_latvec()
            # ------ get structure
            tmp_struc = self._gen_struc_wo_spg()
            if tmp_struc is not None:    # success of generation
                # ------ scale volume
                if self.vol_mu is not None:
                    vol = random.gauss(mu=self.vol_mu, sigma=self.vol_sigma)
                    tmp_struc.scale_lattice(volume=vol)
                    success, mindist_ij, dist = check_distance(tmp_struc,
                                                               self.atype,
                                                               self.mindist)
                    if not success:
                        sys.stderr.write('mindist in gen_wo_spg: {} - {}, {}. retry.\n'.format(
                            self.atype[mindist_ij[0]],
                            self.atype[mindist_ij[1]],
                            dist))
                        continue    # failure
                # --
                # ------ check actual space group using pymatgen
                try:
                    spg_sym, spg_num = tmp_struc.get_space_group_info(
                        symprec=self.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                # ------ register the structure in pymatgen format
                cid = len(init_struc_data) + id_offset
                init_struc_data[cid] = tmp_struc
                print('Structure ID {0:>6} was generated.'
                      ' Space group: {1:>3} {2}'.format(cid, spg_num, spg_sym))
                # ------ save poscar
                if init_pos_path is not None:
                    out_poscar(tmp_struc, cid, init_pos_path)
        self.init_struc_data = init_struc_data

    def gen_with_find_wy(self, nstruc, spgnum='all', id_offset=0,
                         init_pos_path=None, fwpath='./find_wy'):
        '''
        Generate random structures with space gruop information
        using find_wy program

        # ---------- args
        nstruc (int): number of generated structures

        spgnum ('all' or list): space group numbers which you use
                                'all' --> 1--230
                                list --> e.g. [1, 3, 100 ..., 229, 230]

        id_offset (int): structure ID starts from id_offset
                         e.g. nstruc = 3, id_offset = 10
                              you obtain ID 10, ID 11, ID 12

        init_pos_path (str): specify a path of file
                             if you write POSCAR data of init_struc_data
                             ATTENSION: data are appended to the specified file

        fwpath (str): specify a path for a executable file of find_wy program

        # ---------- comment
        generated init_struc_data is saved in self.init_struc_data
        '''
        # ---------- check args
        if type(nstruc) is int and nstruc > 0:
            pass
        else:
            raise ValueError('nstruc must be positive int')
        if spgnum == 'all' or type(spgnum) is list:
            pass
        else:
            raise ValueError('spgnum is wrong. spgnum = {}'.format(spgnum))
        if type(id_offset) is not int:
            raise TypeError('id_offset must be int')
        if init_pos_path is None or type(init_pos_path) is str:
            pass
        else:
            raise ValueError('init_pos_path is wrong.'
                             ' init_pos_path = {}'.format(init_pos_path))
        if not os.path.isfile(fwpath):
            raise IOError('There is no find_wy program in {}'.format(fwpath))
        # ---------- initialize
        init_struc_data = {}
        # ---------- cd tmp_gen_struc
        if not os.path.isdir('tmp_gen_struc'):
            os.mkdir('tmp_gen_struc')
        os.chdir('tmp_gen_struc')

        # ---------- generate structures
        while len(init_struc_data) < nstruc:
            # ------ get spg, a, b, c, alpha, beta, gamma in self.*
            self._gen_lattice(spgnum)
            # ------ get cosa, cosb, and cosg in self.*
            self._calc_cos()
            # ------ write an input file for find_wy
            self._fw_input()
            # ------ loop for same fw_input
            cnt = 0
            while cnt <= self.maxcnt:
                # -- run find_wy
                with open('log_find_wy', 'w') as f:
                    subprocess.call([fwpath, 'input'], stdout=f, stderr=f)
                # -- generate a structure using POS_WY_SKEL_ALL.json
                if not os.path.isfile('POS_WY_SKEL_ALL.json'):
                    wyflag = False
                    break
                wyflag, tmp_struc = self._gen_struc_with_spg()
                if wyflag is False:    # Failure
                    os.remove('POS_WY_SKEL_ALL.json')
                    cnt += 1
                    continue
                else:    # Success
                    self._rm_files()    # rm input POS_WY_SKEL_ALL.json
                    break         # break fw_input loop
            if wyflag is False:
                # -- maximum trial or no POS_WY_SKEL_ALL.json file
                self._rm_files()    # clean
                continue      # to new fw_input
            # ------ scale volume
            if self.vol_mu is not None:
                vol = random.gauss(mu=self.vol_mu, sigma=self.vol_sigma)
                tmp_struc.scale_lattice(volume=vol)
                success, mindist_ij, dist = check_distance(tmp_struc,
                                                           self.atype,
                                                           self.mindist)
                if not success:
                    sys.stderr.write('mindist in gen_with_find_wy: {} - {}, {}. retry.\n'.format(
                        self.atype[mindist_ij[0]],
                        self.atype[mindist_ij[1]],
                        dist))
                    continue    # failure
            # ------ check actual space group using pymatgen
            try:
                spg_sym, spg_num = tmp_struc.get_space_group_info(
                    symprec=self.symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # ------ register the structure in pymatgen format
            cid = len(init_struc_data) + id_offset
            init_struc_data[cid] = tmp_struc
            print('Structure ID {0:>6} was generated.'
                  ' Space group: {1:>3} --> {2:>3} {3}'.format(
                   cid, self.spg, spg_num, spg_sym))
            # ------ save poscar
            if init_pos_path is not None:
                os.chdir('../')    # temporarily go back to ../
                out_poscar(tmp_struc, cid, init_pos_path)
                os.chdir('tmp_gen_struc')
            # ------ clean
            self._rm_files()
        # ---------- go back to ..
        os.chdir('../')
        # ---------- init_struc_data
        self.init_struc_data = init_struc_data

    def _get_atomlist(self):
        '''
        e.g. Na2Cl2
            atomlist = ['Na', 'Na', 'Cl', 'Cl']
        '''
        atomlist = []
        for i in range(len(self.atype)):
            atomlist += [self.atype[i]]*self.nat[i]
        self.atomlist = atomlist

    def _gen_lattice(self, spgnum):
        # ---------- for spgnum = 0: no space group
        if spgnum == 0:
            crystal_systems = ['Triclinic',
                               'Monoclinic',
                               'Orthorhombic',
                               'Tetragonal',
                               'Rhombohedral',
                               'Hexagonal',
                               'Cubic']
            spg = 0
            csys = random.choice(crystal_systems)
        # ---------- for spgnum 1--230
        else:
            # ------ spgnum --> spg
            if spgnum == 'all':
                spg = random.randint(1, 230)
            else:
                spg = random.choice(spgnum)
            if 1 <= spg <= 2:
                csys = 'Triclinic'
            elif 3 <= spg <= 15:
                csys = 'Monoclinic'
            elif 16 <= spg <= 74:
                csys = 'Orthorhombic'
            elif 75 <= spg <= 142:
                csys = 'Tetragonal'
            elif 143 <= spg <= 167:
                # trigonal includes rhombohedral in find_wy
                csys = 'Trigonal'
            elif 168 <= spg <= 194:
                csys = 'Hexagonal'
            elif 195 <= spg <= 230:
                csys = 'Cubic'
            else:
                raise ValueError('spg is wrong')
        # ---------- generate lattice constants a, b, c, alpha, beta, gamma
        if csys == 'Triclinic':
            t1 = random.uniform(self.minlen, self.maxlen)
            t2 = random.uniform(self.minlen, self.maxlen)
            t3 = random.uniform(self.minlen, self.maxlen)
            t = [t1, t2, t3]
            t.sort()
            a, b, c = t
            r = random.random()
            if r < 0.5:    # Type I
                alpha = 90.0 - random.uniform(0, self.dangle)
                beta  = 90.0 - random.uniform(0, self.dangle)
                gamma = 90.0 - random.uniform(0, self.dangle)
            else:    # Type II
                alpha = 90.0 + random.uniform(0, self.dangle)
                beta  = 90.0 + random.uniform(0, self.dangle)
                gamma = 90.0 + random.uniform(0, self.dangle)
        elif csys == 'Monoclinic':
            a = random.uniform(self.minlen, self.maxlen)
            b = random.uniform(self.minlen, self.maxlen)
            c = random.uniform(self.minlen, self.maxlen)
            if a > c:
                a, c = c, a
            alpha = gamma = 90.0
            beta = 90.0 + random.uniform(0, self.dangle)
        elif csys == 'Orthorhombic':
            t1 = random.uniform(self.minlen, self.maxlen)
            t2 = random.uniform(self.minlen, self.maxlen)
            t3 = random.uniform(self.minlen, self.maxlen)
            t = [t1, t2, t3]
            t.sort()
            a, b, c = t
            alpha = beta = gamma = 90.0
        elif csys == 'Tetragonal':
            a = b = random.uniform(self.minlen, self.maxlen)
            c = random.uniform(self.minlen, self.maxlen)
            alpha = beta = gamma = 90.0
        elif csys == 'Trigonal':
            a = b = random.uniform(self.minlen, self.maxlen)
            c = random.uniform(self.minlen, self.maxlen)
            alpha = beta = 90.0
            gamma = 120.0
        elif csys == 'Rhombohedral':
            a = b = c = random.uniform(self.minlen, self.maxlen)
            alpha = beta = gamma = 90 + random.uniform(-self.dangle,
                                                       self.dangle)
        elif csys == 'Hexagonal':
            a = b = random.uniform(self.minlen, self.maxlen)
            c = random.uniform(self.minlen, self.maxlen)
            alpha = beta = 90.0
            gamma = 120.0
        elif csys == 'Cubic':
            a = b = c = random.uniform(self.minlen, self.maxlen)
            alpha = beta = gamma = 90.0
        self.spg = spg
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

    def _calc_latvec(self):
        # ---------- degree to radian
        alpha_rad = math.radians(self.alpha)
        beta_rad = math.radians(self.beta)
        gamma_rad = math.radians(self.gamma)
        # ---------- calculate components
        bx = self.b*math.cos(gamma_rad)
        by = self.b*math.sin(gamma_rad)
        cx = self.c*math.cos(beta_rad)
        cy = (self.c*math.cos(alpha_rad)
              - cx*math.cos(gamma_rad))/math.sin(gamma_rad)
        cz = math.sqrt(self.c*self.c - cx*cx - cy*cy)
        # ---------- lattice vector as list
        self.va = [self.a, 0.0, 0.0]
        self.vb = [bx, by, 0.0]
        self.vc = [cx, cy, cz]

    def _calc_cos(self):
        # ---------- degree to radian
        a_rad = math.radians(self.alpha)
        b_rad = math.radians(self.beta)
        g_rad = math.radians(self.gamma)
        self.cosa = math.cos(a_rad)
        self.cosb = math.cos(b_rad)
        self.cosg = math.cos(g_rad)

    def _gen_struc_wo_spg(self):
        '''
        Success --> return structure data in pymatgen format
        Failure --> return None
        '''
        # ---------- initialize
        cnt = 0
        incoord = []
        # ---------- generate internal coordinates
        while len(incoord) < self.natot:
            tmp_coord = np.random.rand(3)
            incoord.append(tmp_coord)
            tmp_struc = Structure([self.va, self.vb, self.vc],
                                  self.atomlist[:len(incoord)],
                                  incoord)
            success, mindist_ij, dist = check_distance(tmp_struc,
                                                       self.atype,
                                                       self.mindist)
            if not success:
                sys.stderr.write('mindist in _gen_struc_wo_spg: {} - {}, {}. retry.\n'.format(
                    self.atype[mindist_ij[0]],
                    self.atype[mindist_ij[1]],
                    dist))
                incoord.pop(-1)    # cancel
                cnt += 1
                if self.maxcnt < cnt:
                    return None
        return tmp_struc

    def _fw_input(self):
        with open('input', 'w') as f:
            f.write('nspecies {}\n'.format(len(self.atype)))
            f.write('species_name')
            for aa in self.atype:
                f.write('  {}'.format(aa))
            f.write('\n')
            f.write('species_num')
            for i in self.nat:
                f.write('  {}'.format(i))
            f.write('\n')
            f.write('spacegroup  {}\n'.format(self.spg))
            f.write('originchoice  1\n')
            f.write('\n')
            f.write('a  {}\n'.format(self.a))
            f.write('b  {}\n'.format(self.b))
            f.write('c  {}\n'.format(self.c))
            f.write('cosa  {}\n'.format(self.cosa))
            f.write('cosb  {}\n'.format(self.cosb))
            f.write('cosc  {}\n'.format(self.cosg))
            f.write('\n')
            # f.write('selectone true\n')
            f.write('randomseed auto\n')

    def _gen_struc_with_spg(self):
        '''
        Success --> return True, structure data
        Failure --> return False, _
        '''
        # ---------- load POS_WY_SKEL_ALL.json
        with open('POS_WY_SKEL_ALL.json', 'r') as f:
            wydata = json.load(f)
        # ---------- generate structure
        plat = wydata['primitivevector']
        clat = wydata['conventionalvector']
        n_uniq, wydata_eq_atom = self._get_wydata_eq_atom(wydata)
        eq_atomnames = {}
        eq_positions = {}
        # ---------- equivalent atom loop
        for key, value in sorted(n_uniq.items(), key=lambda x: x[1]):
            # ------ distribute eq atoms. first, special (num_uniqvar = 0),
            #            then, others
            cnt = 0
            while True:
                eq_atomnames[key], eq_positions[key] = self._gen_eq_atoms(
                    wydata_eq_atom[key])
                # -- sort in original order
                atomnames = []
                positions = []
                for key_a, value_a in sorted(eq_atomnames.items()):
                    atomnames += eq_atomnames[key_a]
                    positions += eq_positions[key_a]
                # -- Cartesian coordinate; use clat (not plat)
                cart = []
                for p in positions:
                    v = np.zeros(3)
                    for i in range(3):
                        a = np.array(clat[i])
                        v += p[i] * a
                    cart.append(v)
                # -- check minimum distance
                spgstruc = Structure(plat, atomnames, cart,
                                     coords_are_cartesian=True)
                success, mindist_ij, dist = check_distance(spgstruc,
                                                           self.atype,
                                                           self.mindist)
                if not success:
                    sys.stderr.write('mindist in _gen_struc_with_spg: {} - {}, {}. retry.\n'.format(
                        self.atype[mindist_ij[0]],
                        self.atype[mindist_ij[1]],
                        dist))
                    # failure
                    # num_uniqvar = 0 --> value == 0
                    cnt = self.maxcnt + 1 if value == 0 else cnt + 1
                    if self.maxcnt < cnt:
                        return False, spgstruc    # spgstruc is dummy
                else:
                    break    # break while loop --> next eq atoms
        return True, spgstruc

    def _get_wydata_eq_atom(self, wydata):
        i = 0    # count eq_atom, not atom
        n_uniq = {}    # num_uniqvar each eq_atom
        wydata_eq_atom = {}    # wydata each eq_atom
        for specie in wydata['atoms']:
            for wydata2 in specie:    # equivalent atom loop
                n_uniq[i] = wydata2[0]['num_uniqvar']
                wydata_eq_atom[i] = wydata2
                i += 1
        return n_uniq, wydata_eq_atom

    def _gen_eq_atoms(self, wydata2):
        eq_atomnames = []
        eq_positions = []
        rval = np.random.random_sample(3)
        for each in wydata2:
            pos = []
            for ch in each['xyzch']:
                if ch == '-2x':
                    pos.append(-2.0 * rval[0])
                elif ch == '-x+y':
                    pos.append(-rval[0] + rval[1])
                elif ch == '-z':
                    pos.append(-rval[2])
                elif ch == '-y':
                    pos.append(-rval[1])
                elif ch == '-x':
                    pos.append(-rval[0])
                elif ch == '0':
                    pos.append(0.0)
                elif ch == 'x':
                    pos.append(rval[0])
                elif ch == 'y':
                    pos.append(rval[1])
                elif ch == 'z':
                    pos.append(rval[2])
                elif ch == 'x-y':
                    pos.append(rval[0] - rval[1])
                elif ch == '2x':
                    pos.append(2.0 * rval[0])
                else:
                    raise ValueError('unknown ch in conversion in gen_wycoord')
            pos = np.array(pos)
            eq_positions.append(pos + each['add'])
            eq_atomnames.append(each['name'])
        return eq_atomnames, eq_positions

    def _rm_files(self, files=['input', 'POS_WY_SKEL_ALL.json']):
        for rfile in files:
            if os.path.isfile(rfile):
                os.remove(rfile)
