'''
Utility for structures
'''

from logging import getLogger
import os

import numpy as np

from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pyxtal.tolerance import Tol_matrix

from ..IO import read_input as rin


logger = getLogger('cryspy')

def set_mindist(mindist_in, factor, dummy=False, mpi_rank=0):
    # ---------- dummy atom in mol_bs
    if dummy:
        atype = get_atype_dummy()
    else:
        atype = rin.atype

    # ---------- mindist
    if mindist_in is None:
        # ------ Tol matrix
        if rin.struc_mode in ['crystal']:
            tolmat = Tol_matrix(prototype='atomic', factor=factor)
        elif rin.struc_mode in ['mol', 'mol_bs']:
            tolmat = Tol_matrix(prototype='molecular', factor=factor)
        else:
            logger.error('struc_mode is wrong')
            raise SystemExit(1)
        # ------ set mindist
        mindist = []
        for i, itype in enumerate(atype):
            tmp = []
            for j, jtype in enumerate(atype):
                tmp.append(tolmat.get_tol(itype, jtype))
            mindist.append(tmp)
    else:
        tmp_array = factor * np.array(mindist_in)
        mindist = tmp_array.tolist()

    # ---------- log
    if mpi_rank == 0:
        for i, itype in enumerate(atype):
            for j, jtype in enumerate(atype):
                if i <= j:
                    if dummy:
                        logger.info(f'{rin.mol_file[i]} - {rin.mol_file[j]}: {mindist[i][j]}')
                    else:
                        logger.info(f'{itype} - {jtype}: {mindist[i][j]}')
        if rin.struc_mode == 'mol':
            logger.info('When struc_mode is mol (only random structure, not EA part),\n'
            '- tolerance between monoatomic molecules is multiplied by 0.8 inside pyxtal (not printed above)\n'
            '- H-N, H-O, or H-F tolerance is multiplied by 0.9 inside pyxtal (not printed above)')

    # ---------- return
    return mindist


def get_atype_dummy():
    noble_gas = ['Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He']
    if len(rin.nmol) > len(noble_gas):
        logger.error('len(nmol) > len(noble_gas)')
        raise SystemExit(1)
    atype = noble_gas[:len(rin.nmol)]
    return atype


def out_poscar(struc, cid, fpath):
    # ---------- poscar format
    pos = struc.to(fmt='poscar')
    pos = pos.split('\n')
    blank_indx = pos.index('')    # cut unnecessary parts
    pos = pos[:blank_indx]
    pos[0] = 'ID_{}'.format(cid)    # replace with ID
    lines = [line+'\n' for line in pos]

    # ---------- append POSCAR
    with open(fpath, 'a') as f:
        for line in lines:
            f.write(line)


def out_cif(struc, cid, tmp_path, fpath, symprec=0.1):
    # ---------- opt_CIFS
    cif = CifWriter(struc, symprec=symprec)
    cif.write_file(tmp_path+'tmp.cif')

    # ---------- correct title for VESTA
    #                (need to delete '_chemical_formula_sum'. i don't know why)
    with open(tmp_path+'tmp.cif', 'r') as fcif:
        ciflines = fcif.readlines()
    ciflines[1] = 'data_ID_{}\n'.format(cid)
    if ciflines[11][:21] == '_chemical_formula_sum':
        ciflines.pop(11)
    else:
        logger.error('ciflines[11] is not _chemical_formula_sum,'
                         ' have to fix bag')
        raise SystemExit(1)

    # ---------- cif --> opt_cifs
    with open(fpath, 'a') as foptcif:
        for line in ciflines:
            foptcif.write(line)

    # ---------- clean tmp.cif
    os.remove(tmp_path+'tmp.cif')


def frac_coord_zero_one(struc_in):
    '''
    fractional coordinates: 0.0 <= x,y,z < 1.0
    e.g. [0.0, -0.25, 0.7] --> [0.0, 0.75, 0.7]

    # ---------- args
    struc_in: structure data in pymatgen format
    '''
    struc = struc_in.copy()
    for i in range(struc.num_sites):
        struc[i] = struc[i].to_unit_cell()
    return struc


def origin_shift(struc_in):
    '''
    Randomly shift the origin of struc_in

    # ---------- args
    struc_in: structure data in pymatgen format

    # ---------- return
    origin shifted structure (not change original struc_in)
    '''
    struc = struc_in.copy()
    coords_trans = struc.frac_coords + np.random.rand(3)
    struc_shift = Structure(struc.lattice, struc.species, coords_trans)
    struc_shift = frac_coord_zero_one(struc_shift)
    return struc_shift


def sort_by_atype(struc, atype):
    '''
    return a structre sorted by atype order as a new structure
    '''
    return struc.get_sorted_structure(
        key=lambda x: atype.index(x.species_string))


def check_distance(struc, atype, mindist, check_all=False):
    '''
    # ---------- args
    struc: structure data in pymatgen format
    atype (list): e.g. ['Li', 'Co, 'O']
    mindist (2d list) : e.g. [[2.0, 2.0, 1.2],
                              [2.0, 2.0, 1.2],
                              [1.2, 1.2, 1.5]]
    check_all (bool) : if True, check all atom pairs, return dist_list.
                       if False, stop when (dist < mindist) is found,
                                 return True or False (see below)
    # ---------- return
    (check_all=False) True, None, None: nothing smaller than mindist
    (check_all=False) False, (i, j), dist: something smaller than mindst
                                     here, (i, j) means mindist(i, j)
    (check_all=True) dist_list: if dist_list is vacant,
                                    nothing smaller than mindist
    '''

    # ---------- initialize
    if check_all:
        dist_list = []    # [(i, j, dist), (i, j, dist), ...]

    # ---------- in case there is only one atom
    if struc.num_sites == 1:
        dist = min(struc.lattice.abc)
        if dist < mindist[0][0]:
            if check_all:
                dist_list.append((0, 0, dist))
                return dist_list
            return False, (0, 0), dist
        return True, None, None

    # ---------- normal case
    for i in range(struc.num_sites):
        for j in range(i):
            dist = struc.get_distance(j, i)
            i_type = atype.index(struc[i].species_string)
            j_type = atype.index(struc[j].species_string)
            if dist < mindist[i_type][j_type]:
                if check_all:
                    dist_list.append((j, i, dist))
                else:
                    return False, (i_type, j_type), dist

    # ---------- return
    if check_all:
        if dist_list:
            dist_list.sort()    # sort
            return dist_list
        else:
            return dist_list    # dist_list is vacant list
    return True, None, None


def scale_cell_mol(struc, mol_data, volume, dist_digit=3):
    '''
    scale cell without changing molecule size
    for molecular structures generated by PyXtal
    '''
    # ---------- get molecules
    molecules = []
    for mol in mol_data:
        if len(mol.species) > 1:    # mol.species == 1 <-- atom
            molecules.append(mol)

    # ---------- structure of molecule
    mols_cart_coords = []
    mols_frac_Gs = []
    for mol in molecules:
        # ------ start atom for the molecule
        mol_pos = []
        for i, element in enumerate(struc.species):
            # -- check len
            if len(struc.species) < mol.num_sites + i:
                break
            # -- check number of atoms
            count = 0
            for j in range(mol.num_sites):
                if struc.species[i+j] == mol.species[j]:
                    count += 1
            if count == mol.num_sites:
                mol_pos.append(i)
        # ------ calculate the center of gravity
        in_dists = np.linalg.norm(mol.cart_coords - mol.cart_coords[0], axis=1)
        frac_Gs = []
        mol_cart_coords = []
        # -- take periodic boundary into account
        #    here, check distance from the reference atom (mol.cart_coords[0])
        #    to find a correct site for calculating the center of gravity
        for st in mol_pos:
            frac_mol_coords = struc.frac_coords[st:st + mol.num_sites]
            cart_mol_basis = struc.lattice.get_cartesian_coords(frac_mol_coords[0])
            fix_frac_coords = []
            for frac, dist in zip(frac_mol_coords, in_dists):
                all_pattern = [[frac[0]], [frac[1]], [frac[2]]]
                basis_dist = frac - frac_mol_coords[0]
                for axis in range(3):
                    if basis_dist[axis] > 0:
                        all_pattern[axis].append(frac[axis]-1)
                    elif basis_dist[axis] < 0:
                        all_pattern[axis].append(frac[axis]+1)
                # check all patterns
                break_loop = False
                for a_coord in all_pattern[0]:
                    for b_coord in all_pattern[1]:
                        for c_coord in all_pattern[2]:
                            tmp_frac_coord = [a_coord, b_coord, c_coord]
                            tmp_cart_coord = struc.lattice.get_cartesian_coords(tmp_frac_coord)
                            tmp_cart_dist = np.linalg.norm(tmp_cart_coord - cart_mol_basis)
                            if round(dist, dist_digit) == round(tmp_cart_dist, dist_digit):
                                break_loop = True
                                fix_frac_coords.append(tmp_frac_coord)
                                break
                        if break_loop:
                            break
                    if break_loop:
                        break
            if not len(fix_frac_coords) == mol.num_sites:
                return False
            fix_frac_coords = np.array(fix_frac_coords)
            # -- get the center of gravity
            frac_G = fix_frac_coords.sum(axis=0)/mol.num_sites
            frac_Gs.append(frac_G)
            # -- save structure of mol
            fix_cart_coords = struc.lattice.get_cartesian_coords(fix_frac_coords)
            cart_G = fix_cart_coords.sum(axis=0)/mol.num_sites
            mol_cart_coord = fix_cart_coords - cart_G
            mol_cart_coords.append(mol_cart_coord)
        mols_frac_Gs.append(frac_Gs)
        mols_cart_coords.append(mol_cart_coords)
        # ------ remove the molecule from sturc
        for st in reversed(mol_pos):
            struc.remove_sites(range(st, st + mol.num_sites))

    # ---------- scale lattice
    struc.scale_lattice(volume)

    # ---------- put the molecule of original size
    for mol_index in range(len((molecules))):
        molecule = molecules[mol_index]
        for frac_G, cart_mol in zip(mols_frac_Gs[mol_index],
                                    mols_cart_coords[mol_index]):
            cart_G = struc.lattice.get_cartesian_coords(frac_G)
            cart_mol_coord = cart_mol + cart_G
            for spe, spe_coord in zip(molecule.species, cart_mol_coord):
                struc.append(spe, spe_coord, coords_are_cartesian=True)

    # ---------- return new struc
    return struc


def rot_mat(angles, seq='zyx', degree=False):
    '''
    Calculate rotation matrix

    # ---------- input
    angles: in radian (degree=False)
    seq:
    degree:
    '''
    if not len(seq) == len(angles):
        logger.error('not len(seq) == len(angles)')
        raise SystemExit(1)

    if degree:
        angles = np.deg2rad(angles)

    R = np.eye(3)
    angles = angles[::-1]
    seq = seq[::-1]
    for axis, angle in zip(seq, angles):
        if 'x' == axis:
            tmp_R = np.array([[1, 0, 0],
                              [0, np.cos(angle), -np.sin(angle)],
                              [0, np.sin(angle), np.cos(angle)]])
        elif 'y' == axis:
            tmp_R = np.array([[np.cos(angle), 0, np.sin(angle)],
                              [0, 1, 0],
                              [-np.sin(angle), 0, np.cos(angle)]])
        elif 'z' == axis:
            tmp_R = np.array([[np.cos(angle), -np.sin(angle), 0],
                              [np.sin(angle), np.cos(angle), 0],
                              [0, 0, 1]])
        else:
            tmp_R = np.eye(3)
        R = np.matmul(R, tmp_R)

    return R


def find_site(struc, mol_id, group_id, true_dists):
    '''
    find out atom on the unit cell, because don't refer atom in neighbor cell
    '''
    # ---------- init
    fix_frac_coords = []
    fix_group_id = []
    fix_mol_id = []
    grouped_struc_mol_id = []
    id_cnt = 0
    grouped_struc = []
    # ----------  separate molecules per group_id
    for j in range(len(struc)):
        remove_indx = []
        for i, s in enumerate(struc):
            if id_cnt != group_id[i]:
                remove_indx.append(i)
        if len(remove_indx) != len(struc):
            tmp_struc = struc.copy()
            tmp_struc.remove_sites(remove_indx)
            grouped_struc.append(tmp_struc)
            for n in range(len(struc)):
                if group_id[n] == id_cnt:
                    grouped_struc_mol_id.append(mol_id[n])
                    break
        id_cnt += 1
    fix_species = []
    # ---------- compare atom distance to atom distance of mol_data
    for j, ts in enumerate(grouped_struc):
        for i, frac in enumerate(ts.frac_coords):
            dist = true_dists[grouped_struc_mol_id[j]][i]
            all_pattern = [[frac[0]], [frac[1]], [frac[2]]]
            basis_dist = frac - ts[0].frac_coords
            for axis in range(3):
                if basis_dist[axis] >= 0:
                    all_pattern[axis].append(frac[axis]-1)
                elif basis_dist[axis] < 0:
                    all_pattern[axis].append(frac[axis]+1)
            # ------ check all patterns
            # break_loop = False
            min_dist = 1000000
            for a_coord in all_pattern[0]:
                for b_coord in all_pattern[1]:
                    for c_coord in all_pattern[2]:
                        tmp_frac_coord = [a_coord, b_coord, c_coord]
                        tmp_cart_coord = struc.lattice.get_cartesian_coords(tmp_frac_coord)
                        base_cart = struc.lattice.get_cartesian_coords(ts.frac_coords[0])
                        tmp_cart_dist = np.linalg.norm(tmp_cart_coord - base_cart)
                        dist_rmse = np.sqrt(np.mean((dist - tmp_cart_dist)**2))
                        if dist_rmse < min_dist:
                            min_dist = dist_rmse
                            min_dist_frac_coord = tmp_frac_coord
            fix_frac_coords.append(min_dist_frac_coord)
            fix_group_id.append(j)
            for n, gid in enumerate(group_id):
                if gid == j:
                    fix_mol_id.append(grouped_struc_mol_id[j])
                    break
            fix_species.append(grouped_struc[j][i].specie)
    return np.array(fix_frac_coords), fix_group_id, fix_species, fix_mol_id


def cal_g(struc, mol_id, group_id, true_dists):
    # ----------  calculate molecule center
    mol_g_A = []
    frac_coords, fix_group_id, _, _ = find_site(struc, mol_id, group_id, true_dists)
    for cnt in range(len(struc)):
        acnt = 0
        cal_g = 0
        for i, fc in enumerate(frac_coords):
            if cnt == fix_group_id[i]:
                acnt += 1
                cal_g += fc
        if acnt == 0:
            break
        mol_g_A.append(cal_g/acnt)
    return mol_g_A


def sort_by_atype_mol(struc, atype, mol_id, group_id):
    # ---------- work the same as sort_by_atype, and sort mol_id and group_id
    copy_struc = struc.copy()
    sorted_struc = sort_by_atype(struc, atype)
    sorted_mol_id = []
    sorted_group_id = []
    for i in range(len(sorted_struc)):
        for j in range(len(copy_struc)):
            if np.allclose(sorted_struc.frac_coords[i], copy_struc.frac_coords[j]):
                sorted_mol_id.append(mol_id[j])
                sorted_group_id.append(group_id[j])
    return sorted_struc, sorted_mol_id, sorted_group_id


def get_nat(struc, atype):
    compos = struc.composition.as_dict()
    nat = [int(compos[at]) for at in atype]
    ratio = [x/sum(nat) for x in nat]
    return nat, ratio