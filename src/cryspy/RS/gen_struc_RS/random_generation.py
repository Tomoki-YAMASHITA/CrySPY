'''
Random structure generation w/o pyxtal
'''

import json
from logging import getLogger
import math
import os
import random
import subprocess

import numpy as np
from pymatgen.core import Structure

from ...util.struc_util import check_distance


logger = getLogger('cryspy')


def gen_wo_spg(rin, nstruc, mindist, id_offset=0, vc=False):
    '''
    Generate random structures without space group information

    # ---------- args
    rin (instance of ReadInput): input parameters
    nstruc (int): number of structures to be generated
    mindist (): minimum interatomic distance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                            you obtain ID 10, ID 11, ID 12
    vc (bool): variable composition. it needs ll_nat and ul_nat

    # ---------- return
    init_struc_data (dict): {ID: pymatgen structure data}
    '''

    # ---------- initialize
    init_struc_data = {}

    # ---------- generate structures
    while len(init_struc_data) < nstruc:
        # ------ vc
        if not vc:
            numIons = rin.nat
        else:    # variable composition
            numIons = [random.randint(l, u) for l, u in zip(rin.ll_nat, rin.ul_nat)]
            numIons = tuple(numIons)
        atomlist = _get_atomlist(rin, numIons)
        # ------ get spg, a, b, c, alpha, beta, gamma
        spg, a, b, c, alpha, beta, gamma = _gen_lattice(rin)
        # ------ get a1, a2, a3
        a1, a2, a3 = _calc_latvec(a, b, c, alpha, beta, gamma)
        # ------ get structure
        tmp_struc = _gen_struc_wo_spg(rin, numIons, atomlist, a1, a2, a3, mindist)
        if tmp_struc is not None:    # success of generation
            # ------ scale volume
            if rin.vol_mu is not None:
                vol = random.gauss(mu=rin.vol_mu, sigma=rin.vol_sigma)
                tmp_struc.scale_lattice(volume=vol)
                success, mindist_ij, dist = check_distance(tmp_struc,
                                                            rin.atype,
                                                            mindist)
                if not success:
                    type0 = rin.atype[mindist_ij[0]]
                    type1 = rin.atype[mindist_ij[1]]
                    logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
                    continue    # failure
            # ------ check actual space group using pymatgen
            try:
                spg_sym, spg_num = tmp_struc.get_space_group_info(
                    symprec=rin.symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # ------ register the structure in pymatgen format
            cid = len(init_struc_data) + id_offset
            init_struc_data[cid] = tmp_struc
            logger.info(f'Structure ID {cid:>6}: {numIons}'
                    f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return init_struc_data


def gen_with_find_wy(rin, nstruc, mindist, id_offset=0,
                        fwpath='find_wy', mpi_rank=0, vc=False):
    '''
    Generate random structures with space gruop information
    using find_wy program

    # ---------- args
    rin (instance of ReadInput): input parameters
    nstruc (int): number of generated structures
    mindist (): minimum interatomic distance
    id_offset (int): structure ID starts from id_offset
                        e.g. nstruc = 3, id_offset = 10
                            you obtain ID 10, ID 11, ID 12
    fwpath (str): specify a path for a executable file of find_wy program
    mpi_rank (int): rank of MPI process
    vc (bool): variable composition. it needs ll_nat and ul_nat

    # ---------- return
    init_struc_data (dict): {ID: pymatgen Structre}
    '''

    # ---------- initialize
    init_struc_data = {}

    # ---------- cd tmp_gen_struc
    os.makedirs(f'tmp_gen_struc/rank_{mpi_rank}', exist_ok=True)
    os.chdir(f'tmp_gen_struc/rank_{mpi_rank}')

    # ---------- generate structures
    while len(init_struc_data) < nstruc:
        # ------ vc
        if not vc:
            numIons = rin.nat
        else:    # variable composition
            numIons = [random.randint(l, u) for l, u in zip(rin.ll_nat, rin.ul_nat)]
            numIons = tuple(numIons)
        # ------ get spg, a, b, c, alpha, beta, gamma
        spg, a, b, c, alpha, beta, gamma = _gen_lattice(rin)
        # ------ get cosa, cosb, and cosc
        cosa, cosb, cosg = _calc_cos(alpha, beta, gamma)
        # ------ write an input file for find_wy
        _fw_input(rin, numIons, spg, a, b, c, cosa, cosb, cosg)
        # ------ loop for same fw_input
        cnt = 0
        while cnt <= rin.maxcnt:
            # -- run find_wy
            with open('log_find_wy', 'w') as f:
                subprocess.call([fwpath, 'input'], stdout=f, stderr=f)
            # -- generate a structure using POS_WY_SKEL_ALL.json
            if not os.path.isfile('POS_WY_SKEL_ALL.json'):
                wyflag = False
                break
            wyflag, tmp_struc = _gen_struc_with_spg(rin, mindist)
            if wyflag is False:    # Failure
                os.remove('POS_WY_SKEL_ALL.json')
                cnt += 1
                continue
            else:    # Success
                _rm_files()    # rm input POS_WY_SKEL_ALL.json
                break          # break fw_input loop
        if wyflag is False:
            # -- maximum trial or no POS_WY_SKEL_ALL.json file
            _rm_files()    # clean
            continue       # to new fw_input
        # ------ scale volume
        if rin.vol_mu is not None:
            vol = random.gauss(mu=rin.vol_mu, sigma=rin.vol_sigma)
            tmp_struc.scale_lattice(volume=vol)
            success, mindist_ij, dist = check_distance(tmp_struc,
                                                        rin.atype,
                                                        mindist)
            if not success:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
                continue    # failure
        # ------ check actual space group using pymatgen
        try:
            spg_sym, spg_num = tmp_struc.get_space_group_info(
                symprec=rin.symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ------ tmp_struc --> init_struc_data
        cid = len(init_struc_data) + id_offset
        init_struc_data[cid] = tmp_struc
        logger.info(f'Structure ID {cid:>6}: {numIons}'
                f' Space group: {spg:>3} --> {spg_num:>3} {spg_sym}')
        # ------ clean
        _rm_files()

    # ---------- go back to ..
    os.chdir('../../')

    # ---------- return
    return init_struc_data


def _get_atomlist(rin, numIons):
    '''
    e.g. Na2Cl2
        atomlist = ['Na', 'Na', 'Cl', 'Cl']
    '''
    atomlist = []
    for i in range(len(rin.atype)):
        atomlist += [rin.atype[i]]*numIons[i]
    return atomlist

def _gen_lattice(rin):
    # ---------- for spgnum = 0: no space group
    if rin.spgnum == 0:
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
        if rin.spgnum == 'all':
            spg = random.randint(1, 230)
        else:
            spg = random.choice(rin.spgnum)
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
            logger.error('spg is wrong')
            raise SystemExit(1)
    # ---------- generate lattice constants a, b, c, alpha, beta, gamma
    if csys == 'Triclinic':
        t1 = random.uniform(rin.minlen, rin.maxlen)
        t2 = random.uniform(rin.minlen, rin.maxlen)
        t3 = random.uniform(rin.minlen, rin.maxlen)
        t = [t1, t2, t3]
        t.sort()
        a, b, c = t
        r = random.random()
        if r < 0.5:    # Type I
            alpha = 90.0 - random.uniform(0, rin.dangle)
            beta  = 90.0 - random.uniform(0, rin.dangle)
            gamma = 90.0 - random.uniform(0, rin.dangle)
        else:    # Type II
            alpha = 90.0 + random.uniform(0, rin.dangle)
            beta  = 90.0 + random.uniform(0, rin.dangle)
            gamma = 90.0 + random.uniform(0, rin.dangle)
    elif csys == 'Monoclinic':
        a = random.uniform(rin.minlen, rin.maxlen)
        b = random.uniform(rin.minlen, rin.maxlen)
        c = random.uniform(rin.minlen, rin.maxlen)
        if a > c:
            a, c = c, a
        alpha = gamma = 90.0
        beta = 90.0 + random.uniform(0, rin.dangle)
    elif csys == 'Orthorhombic':
        t1 = random.uniform(rin.minlen, rin.maxlen)
        t2 = random.uniform(rin.minlen, rin.maxlen)
        t3 = random.uniform(rin.minlen, rin.maxlen)
        t = [t1, t2, t3]
        t.sort()
        a, b, c = t
        alpha = beta = gamma = 90.0
    elif csys == 'Tetragonal':
        a = b = random.uniform(rin.minlen, rin.maxlen)
        c = random.uniform(rin.minlen, rin.maxlen)
        alpha = beta = gamma = 90.0
    elif csys == 'Trigonal':
        a = b = random.uniform(rin.minlen, rin.maxlen)
        c = random.uniform(rin.minlen, rin.maxlen)
        alpha = beta = 90.0
        gamma = 120.0
    elif csys == 'Rhombohedral':
        a = b = c = random.uniform(rin.minlen, rin.maxlen)
        alpha = beta = gamma = 90 + random.uniform(-rin.dangle,
                                                    rin.dangle)
    elif csys == 'Hexagonal':
        a = b = random.uniform(rin.minlen, rin.maxlen)
        c = random.uniform(rin.minlen, rin.maxlen)
        alpha = beta = 90.0
        gamma = 120.0
    elif csys == 'Cubic':
        a = b = c = random.uniform(rin.minlen, rin.maxlen)
        alpha = beta = gamma = 90.0

    # ---------- return
    return spg, a, b, c, alpha, beta, gamma


def _calc_latvec(a, b, c, alpha, beta, gamma):
    # ---------- degree to radian
    alpha_rad = math.radians(alpha)
    beta_rad = math.radians(beta)
    gamma_rad = math.radians(gamma)
    # ---------- calculate components
    bx = b*math.cos(gamma_rad)
    by = b*math.sin(gamma_rad)
    cx = c*math.cos(beta_rad)
    cy = (c*math.cos(alpha_rad)
            - cx*math.cos(gamma_rad))/math.sin(gamma_rad)
    cz = math.sqrt(c*c - cx*cx - cy*cy)
    # ---------- lattice vector as list
    a1 = [a, 0.0, 0.0]
    a2 = [bx, by, 0.0]
    a3 = [cx, cy, cz]
    # ---------- return
    return a1, a2, a3


def _calc_cos(alpha, beta, gamma):
    # ---------- degree to radian
    a_rad = math.radians(alpha)
    b_rad = math.radians(beta)
    g_rad = math.radians(gamma)
    cosa = math.cos(a_rad)
    cosb = math.cos(b_rad)
    cosg = math.cos(g_rad)
    # ---------- return
    return cosa, cosb, cosg


def _gen_struc_wo_spg(rin, numIons, atomlist, a1, a2, a3, mindist):
    '''
    Success --> return structure data in pymatgen format
    Failure --> return None
    '''
    # ---------- initialize
    cnt = 0
    incoord = []
    natot = sum(numIons)    # do not use rin.natot because of vc
    # ---------- generate internal coordinates
    while len(incoord) < natot:
        tmp_coord = np.random.rand(3)
        incoord.append(tmp_coord)
        tmp_struc = Structure([a1, a2, a3],
                                atomlist[:len(incoord)],
                                incoord)
        success, mindist_ij, dist = check_distance(tmp_struc,
                                                    rin.atype,
                                                    mindist)
        if not success:
            type0 = rin.atype[mindist_ij[0]]
            type1 = rin.atype[mindist_ij[1]]
            logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
            incoord.pop(-1)    # cancel
            cnt += 1
            if rin.maxcnt < cnt:
                return None
    return tmp_struc


def _fw_input(rin, numIons, spg, a, b, c, cosa, cosb, cosg):
    with open('input', 'w') as f:
        f.write(f'nspecies {len(rin.atype)}\n')
        f.write('species_name')
        for aa in rin.atype:
            f.write(f'  {aa}')
        f.write('\n')
        f.write('species_num')
        for i in numIons:
            f.write(f'  {i}')
        f.write('\n')
        f.write(f'spacegroup  {spg}\n')
        f.write('originchoice  1\n')
        f.write('\n')
        f.write(f'a  {a}\n')
        f.write(f'b  {b}\n')
        f.write(f'c  {c}\n')
        f.write(f'cosa  {cosa}\n')
        f.write(f'cosb  {cosb}\n')
        f.write(f'cosc  {cosg}\n')
        f.write('\n')
        # f.write('selectone true\n')
        f.write('randomseed auto\n')


def _gen_struc_with_spg(rin, mindist):
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
    n_uniq, wydata_eq_atom = _get_wydata_eq_atom(wydata)
    eq_atomnames = {}
    eq_positions = {}
    # ---------- equivalent atom loop
    for key, value in sorted(n_uniq.items(), key=lambda x: x[1]):
        # ------ distribute eq atoms. first, special (num_uniqvar = 0),
        #            then, others
        cnt = 0
        while True:
            eq_atomnames[key], eq_positions[key] = _gen_eq_atoms(
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
                                                        rin.atype,
                                                        mindist)
            if not success:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(f'mindist: {type0} - {type1}, {dist}. retry.')
                # failure
                # num_uniqvar = 0 --> value == 0
                cnt = rin.maxcnt + 1 if value == 0 else cnt + 1
                if rin.maxcnt < cnt:
                    return False, spgstruc    # spgstruc is dummy
            else:
                break    # break while loop --> next eq atoms
    return True, spgstruc

def _get_wydata_eq_atom(wydata):
    i = 0    # count eq_atom, not atom
    n_uniq = {}    # num_uniqvar each eq_atom
    wydata_eq_atom = {}    # wydata each eq_atom
    for specie in wydata['atoms']:
        for wydata2 in specie:    # equivalent atom loop
            n_uniq[i] = wydata2[0]['num_uniqvar']
            wydata_eq_atom[i] = wydata2
            i += 1
    return n_uniq, wydata_eq_atom

def _gen_eq_atoms(wydata2):
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
                logger.error('unknown ch in conversion in gen_wycoord')
                raise SystemExit(1)
        pos = np.array(pos)
        eq_positions.append(pos + each['add'])
        eq_atomnames.append(each['name'])
    return eq_atomnames, eq_positions

def _rm_files(files=['input', 'POS_WY_SKEL_ALL.json']):
    for rfile in files:
        if os.path.isfile(rfile):
            os.remove(rfile)
