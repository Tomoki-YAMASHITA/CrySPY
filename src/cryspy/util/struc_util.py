'''
Utility for structures
'''

from logging import getLogger
import os

import numpy as np
from pymatgen.core import Structure, Molecule
from pyxtal.database.collection import Collection
from pyxtal.tolerance import Tol_matrix

# ---------- import later
#from ..IO.pkl_data import load_cn_comb_data, save_cn_comb_data


logger = getLogger('cryspy')


def set_mindist(atype, mindist_in, factor, struc_mode='crystal',
                dummy=False, mol_file=None, mpi_rank=0, no_print=False):
    # ---------- dummy atom in mol_bs
    if dummy:
        atype = get_atype_dummy(len(mol_file))

    # ---------- mindist
    if mindist_in is None:
        # ------ Tol matrix
        if struc_mode in ['crystal']:
            tolmat = Tol_matrix(prototype='atomic', factor=factor)
        elif struc_mode in ['mol', 'mol_bs']:
            tolmat = Tol_matrix(prototype='molecular', factor=factor)
        else:
            logger.error('struc_mode is wrong')
            try:
                os.remove('lock_cryspy')
            except FileNotFoundError:
                pass
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
                    if not no_print:
                        if dummy:
                            logger.info(f'{mol_file[i]} - {mol_file[j]}: {mindist[i][j]}')
                        else:
                            logger.info(f'{itype} - {jtype}: {mindist[i][j]}')
        if struc_mode == 'mol':
            if not no_print:
                logger.info('When struc_mode is mol (only random structure, not EA part),\n'
                '- tolerance between monoatomic molecules is multiplied by 0.8 inside pyxtal (not printed above)\n'
                '- H-N, H-O, or H-F tolerance is multiplied by 0.9 inside pyxtal (not printed above)')

    # ---------- return
    return mindist


def get_atype_dummy(n_mol_file):
    noble_gas = ['Rn', 'Xe', 'Kr', 'Ar', 'Ne', 'He']
    if n_mol_file > len(noble_gas):
        logger.error('len(mol_file) > len(noble_gas)')
        try:
            os.remove('lock_cryspy')
        except FileNotFoundError:
            pass
        raise SystemExit(1)
    atype = noble_gas[:n_mol_file]
    return atype


def get_mol_data(mol_file):
    '''
    get molecular data

    Input:
        mol_file (tuple)
            e.g. ('Li.xyz', 'PS4.xyz')

    return:
        mol_data (list)
            Molecule data in pymatgen format
    '''
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
            logger.error('no molecular files')
            try:
                os.remove('lock_cryspy')
            except FileNotFoundError:
                pass
            raise SystemExit(1)
        mol_data.append(mol)
    return mol_data


def out_poscar(struc_data: dict, fpath: str, mode='a'):
    '''
    # ---------- args
    struc_data (dict): {'ID': Structure, 'ID': Structure, ...}
    fpath (str): file path for output
    mode (str): 'w' or 'a'
    '''
    # ---------- if mode='w', clear file
    with open(fpath, mode):
        pass
    # ---------- wrtie POSCAR
    for cid, struc in struc_data.items():
        # ---------- poscar format
        pos = struc.to(fmt='poscar')
        pos = pos.split('\n')
        blank_indx = pos.index('')    # cut unnecessary parts
        pos = pos[:blank_indx]
        pos[0] = f'ID_{cid}'    # replace with ID
        lines = [line+'\n' for line in pos]
        # ---------- always mode='a'
        with open(fpath, 'a') as f:
            for line in lines:
                f.write(line)


def out_cif(struc, cid, fpath, symprec=0.01):
    # ---------- str
    str_struc = struc.to(fmt='cif', symprec=symprec)
    # ---------- replace for title in VESTA
    lines = str_struc.split('\n')
    for i, line in enumerate(lines):
        if line.startswith('_chemical_formula_sum'):
            lines[i] = f"_chemical_formula_sum   'ID_{cid}'"
            break
    str_struc = '\n'.join(lines)
    # ---------- cif --> opt_cifs
    with open(fpath, 'a') as foptcif:
        foptcif.write(str_struc)


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


def origin_shift(struc_in, rng=None):
    '''
    Randomly shift the origin of struc_in

    # ---------- args
    struc_in: structure data in pymatgen format

    # ---------- return
    origin shifted structure (not change original struc_in)
    '''
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()
    # ---------- shift
    struc = struc_in.copy()
    coords_trans = struc.frac_coords + rng.random(3)
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
    # ---------- description
    check distance between atoms in struc
    Even if the structure falls back to a lower-dimensional system
    (e.g., from a ternary system with atype = ['Li', 'Co', 'O'] to a binary Li-O system),
    the code can still handle it by checking the atomic species present in the structure.
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
        try:
            os.remove('lock_cryspy')
        except FileNotFoundError:
            pass
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
    nat = tuple([int(compos.get(at, 0)) for at in atype])
    return nat


def remove_zero(atype_in, nat_in, mindist_in):
    '''
    nat_in = (4, 0, 2) --> nat_out = (4, 2)
    atype_in = ('Li', 'Co', 'O') --> atype_out = ('Li', 'O')
    mindist_in = [[1.0, 2.0, 3.0], --> mindist_out = [[1.0, 3.0],
                  [2.0, 1.5, 1.8],                    [3.0, 1.9]]
                  [3.0, 1.8, 1.9]]
    '''

    # ---------- remove from nat
    nat_out = tuple([value for value in nat_in if value != 0])

    # ---------- remove from atype
    atype_out = tuple(
        element for i, element in enumerate(atype_in) if nat_in[i] != 0
    )

    # ---------- remove from mindist
    mindist_out = tuple([
        [value for j, value in enumerate(row) if nat_in[j] != 0]
        for i, row in enumerate(mindist_in) if nat_in[i] != 0
    ])

    # ---------- return
    return atype_out, nat_out, mindist_out


#
# ---------- charge neutrality
#

def calc_cn_comb(ll_nat, ul_nat, charge, use_pkl=True):
    """
    Calculate charge-neutral combinations of atoms

    # ---------- args
    ll_nat (tuple): lower limit of the number of atoms
    ul_nat (tuple): upper limit of the number of atoms
    charge (tuple): charge of each atom type
    use_pkl (bool): if True, load/save data from/to pkl file

    # ---------- retrun
    cn_comb (np.ndarray): array of charge-neutral combinations

    # ---------- example
    e.g. Na-Cl
    ll_nat = (0, 0)
    ul_nat = (4, 4)
    charge = (1, -1)

    cn_comb = array([[1, 1],
                     [2, 2],
                     [3, 3],
                     [4, 4]])

    # ---------- multi-valence example
    e.g. Fe-O
    ll_nat = (4, 4)
    ul_nat = (10, 10)
    charge = ((2, 3), -2)   # (Fe2+, Fe3+), O2-

    Expanded representation:
    ll_ext = (0, 0, 4)      # multivalent: no ll per valence; ll applies to total only
    ul_ext = (10, 10, 10)
    charge_ext = (2, 3, -2)
    """
    # ---------- try loading pkl data if it exists
    if use_pkl:
        from ..IO.pkl_data import load_cn_comb_data, save_cn_comb_data
        try:
            old_ll, old_ul, old_charge, old_cn_comb = load_cn_comb_data()
            if (ll_nat == old_ll) and (ul_nat == old_ul) and (charge == old_charge):
                logger.info('Using old charge neutral combinations data')
                return old_cn_comb
        except Exception:
            pass

    # ---------- convert bounds to np.array
    ll = np.array(ll_nat, dtype=int)
    ul = np.array(ul_nat, dtype=int)
    k_orig = len(charge)

    # ---------- case 1: all charges are single-valence (int)
    if all(isinstance(c, int) for c in charge):
        charges = np.array(charge, dtype=int)
        cn_comb = _calc_cn_comb_single(ll, ul, charges)

    # ---------- case 2: some species have multiple valences
    # expand each valence to a "pseudo-species"
    else:
        valence_list = []        # charges in extended space
        orig_index = []    # mapping: extended index -> original species index
        ll_ext_list = []     # lower bounds in extended space
        for i, c in enumerate(charge):
            if isinstance(c, int):
                # single-valence species: lower bound can be kept
                valence_list.append(c)
                orig_index.append(i)
                ll_ext_list.append(ll[i])
            else:
                # multivalent: each valence may be absent
                for v in c:
                    valence_list.append(v)
                    orig_index.append(i)
                    ll_ext_list.append(0)   # no ll per valence; ll applies to total only
        # ------ convert to np.array
        charges_ext = np.array(valence_list, dtype=int)
        ll_ext = np.array(ll_ext_list, dtype=int)
        ul_ext = np.array([ul[i] for i in orig_index], dtype=int)

        # ---------- charge-neutral combinations in the extended space
        cn_ext = _calc_cn_comb_single(ll_ext, ul_ext, charges_ext)
        if cn_ext.size == 0:
            return np.empty((0, k_orig), dtype=int)

        # ---------- aggregate back to original species
        totals = np.zeros((cn_ext.shape[0], k_orig), dtype=int)    # initialized
        for j, p in enumerate(orig_index):
            totals[:, p] += cn_ext[:, j]

        # ---------- apply original bounds ll_nat, ul_nat
        mask_ll = np.all(totals >= ll, axis=1)
        mask_ul = np.all(totals <= ul, axis=1)
        final_mask = mask_ll & mask_ul
        if not np.any(final_mask):
            return np.empty((0, k_orig), dtype=int)
        totals_valid = totals[final_mask]

        # ---------- remove duplicate compositions
        cn_comb = np.unique(totals_valid, axis=0)

    # ---------- save
    if use_pkl:
        cn_comb_data = (ll_nat, ul_nat, charge, cn_comb)
        save_cn_comb_data(cn_comb_data)
        logger.info(f'Charge neutral combinations saved: {len(cn_comb)}')

    # ---------- return charge neutral combinations
    return cn_comb


def _calc_cn_comb_single(ll, ul, charges):
    """
    Core routine for charge-neutral combinations with single valence per species.

    Parameters
    ----------
    ll : np.ndarray, shape (k,)
        Lower limits of the number of atoms for each species (int).
    ul : np.ndarray, shape (k,)
        Upper limits of the number of atoms for each species (int).
    charges : np.ndarray, shape (k,)
        Charges of each species (single valence only, all ints).

    Returns
    -------
    cn_comb : np.ndarray, shape (N, k)
        Array of charge-neutral combinations.
        Each row is (n_1, ..., n_k).
    """

    # ---------- size
    k = len(charges)

    # ---------- choose one "dependent" species whose charge != 0
    #            Select the dependent species: charge with largest absolute value
    nonzero_q_idx = np.where(charges != 0)[0]
    dep_idx = nonzero_q_idx[np.argmax(np.abs(charges[nonzero_q_idx]))]    # dependent index
    indep_idx = [i for i in range(k) if i != dep_idx]
    q_dep = charges[dep_idx]
    q_indep = charges[indep_idx]

    # ---------- generate combinations for independent variables only
    ranges = [np.arange(ll[i], ul[i] + 1) for i in indep_idx]
    mesh = np.meshgrid(*ranges, indexing='ij')
    comb_indep = np.stack(mesh, axis=-1).reshape(-1, len(indep_idx))  # (M, k-1)

    # ---------- compute dependent variable from charge neutrality
    # q_dep * n_dep + sum(q_indep * n_indep) = 0
    # → n_dep = - sum(q_indep * n_indep) / q_dep
    sum_indep = comb_indep @ q_indep  # (M,)
    divisible_mask = (-sum_indep % q_dep == 0)    # divisible by q_dep
    n_dep = -sum_indep // q_dep

    # ---------- mask
    in_range_mask = (n_dep >= ll[dep_idx]) & (n_dep <= ul[dep_idx])    # ll <= n_dep <= ul
    total_atoms = comb_indep.sum(axis=1) + n_dep    # \sum n_indep + n_dep
    nonzero_mask = (total_atoms != 0)    # total atoms != 0
    final_mask = divisible_mask & in_range_mask & nonzero_mask

    # ---------- construct final combinations
    if not np.any(final_mask):
        cn_comb = np.empty((0, k), dtype=int)    # no valid combinations --> empty array
    else:
        valid_indep = comb_indep[final_mask]
        valid_dep = n_dep[final_mask]
        # ------ reconstruct full combinations
        cn_comb = np.zeros((len(valid_dep), k), dtype=int)    # vacant array
        cn_comb[:, dep_idx] = valid_dep
        for j, i in enumerate(indep_idx):
            cn_comb[:, i] = valid_indep[:, j]

    # ---------- return
    return cn_comb


#
# ---------- composition constraints
#
def precompute_feasible_N(ll_nat, ul_nat, min_comp, max_comp):
    """
    Precompute feasible total atom numbers N under composition constraints.

    Parameters
    ----------
    ll_nat : tuple[int]
        Lower bounds of atom counts for each species.
    ul_nat : tuple[int]
        Upper bounds of atom counts for each species.
    min_comp : tuple[float]
        Lower bounds of composition fractions x_i (for each species).
    max_comp : tuple[float]
        Upper bounds of composition fractions x_i (for each species).

    Returns
    -------
    feasible_N : list of (N, lower, upper)
        Each entry corresponds to one feasible total atom number N.
        - N      : int, total number of atoms
        - lower  : np.ndarray, per-species lower bounds n_i for this N
        - upper  : np.ndarray, per-species upper bounds n_i for this N
    """
    # ---------- convert to np.array
    ll = np.array(ll_nat, dtype=int)
    ul = np.array(ul_nat, dtype=int)
    min_comp = np.array(min_comp, dtype=float)
    max_comp = np.array(max_comp, dtype=float)

    # ---------- lower and upper limits for total atom number N
    N_l = int(ll.sum())
    N_u = int(ul.sum())

    # ---------- find feasible N
    feasible_N = []
    for N in range(max(1, N_l), N_u + 1):    # N must be at least 1
        lower = np.maximum(ll, np.ceil(min_comp * N).astype(int))    # L_i(N)
        upper = np.minimum(ul, np.floor(max_comp * N).astype(int))   # U_i(N)
        # ------ feasibility check: lower_i <= upper_i for all i
        if (lower > upper).any():
            continue
        # ------ existence condition: sum(lower) <= N <= sum(upper)
        if lower.sum() <= N <= upper.sum():
            feasible_N.append((N, lower, upper))

    # ---------- return
    return feasible_N


def sample_nat_from_feasible_N(feasible_N, rng=None):
    """
    Parameters
    ----------
    feasible_N : list[(N, lower, upper)] or None
        Output of precompute_feasible_N. If None, it will be computed inside.
    rng : np.random.Generator or None
        Random number generator. If None, np.random.default_rng() is used.

    Returns
    -------
    nat : np.ndarray, shape (k,)
        One sampled integer composition (n_1, ..., n_k).
    N : int
        Total atom number used for this sample.
    """
    # ---------- rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- choose one (N, lower, upper) at random
    N, lower, upper = feasible_N[rng.integers(len(feasible_N))]

    # ---------- sample n_i within [lower_i, upper_i] with sum_i n_i = N
    nat = lower.copy()    # initialize with lower bounds
    remaining = int(N - nat.sum())
    capacity = (upper - lower).astype(int)
    k = len(nat)
    # ------ suffix capacity
    # suffix_capacity[i] = sum_{j=i}^{k-1} capacity[j]
    suffix_capacity = np.zeros(k + 1, dtype=int)
    for i in range(k - 1, -1, -1):
        suffix_capacity[i] = suffix_capacity[i + 1] + capacity[i]
    # ------ sample one by one
    for i in range(k - 1):
        max_add = min(capacity[i], remaining)
        min_add = max(0, remaining - suffix_capacity[i + 1])
        add_i = rng.integers(min_add, max_add + 1)
        nat[i] += add_i
        remaining -= add_i
    # ------ last species takes the rest (guaranteed to fit)
    nat[-1] += remaining

    # ---------- return
    return nat, N