from logging import getLogger

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie

from ...util.struc_util import sort_by_atype, get_nat
from ...util.struc_util import sort_by_atype_mol, check_distance, cal_g, find_site


logger = getLogger('cryspy')


def gen_strain(
        atype,
        mindist,
        struc_data,
        sp,
        n_strain,
        id_start=None,
        symprec=0.01,
        sigma_st=0.5,
        maxcnt_ea=50,
        struc_mol_id=None,
        molecular=False,
        protect_mol_struc=True,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_strain (int): number of structures to generate by strain
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    sigma_st (float): standard deviation for strain matrix
    maxcnt_ea (int): maximum number of trial in permutation

    # ---------- return
    children (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'strain'}
    '''
    # ---------- initialize
    struc_cnt = 0
    children = {}
    #children_mol_id = {}
    parents = {}
    operation = {}

    # ---------- id_offset
    if id_start is None:
        cid = max(struc_data.keys()) + 1
    else:
        if id_start < (max(struc_data.keys()) + 1):
            logger.error('id_start is already included in structure ID of the data')
        else:
            cid = id_start

    # ---------- generate structures by strain
    while struc_cnt < n_strain:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ generate child
        if molecular:
            logger.error('protect_mol_struc is not implemented yet.')
            # if protect_mol_struc:
                # child, mol_id = gen_child_mol(rin,
                #                                     struc_data[pid],
                #                                     struc_mol_id[pid])
            #else:
            #    child = gen_child(atype, mindist, parent_struc, sigma_st, maxcnt_ea)
        else:
            child = gen_child(atype, mindist, parent_A, sigma_st, maxcnt_ea)
        # ------ success
        if child is not None:
            children[cid] = child
            # if molecular:
            #     children_mol_id[cid] = mol_id
            parents[cid] = (pid_A, )    # tuple
            operation[cid] = 'strain'
            try:
                spg_sym, spg_num = child.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(child, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                    f' from {pid_A:>6} by strain.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return children, parents, operation


def gen_child(atype, mindist, parent_A, sigma_st=0.5, maxcnt_ea=50):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    sigma_st (float): standard deviation for strain matrix
    maxcnt_ea (int): maximum number of trial in crossover

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''
    # ---------- initialize
    child = parent_A.copy()    # keep original structure
    lat_mat = child.lattice.matrix.T    # lattice vector as matrix
    cnt = 0

    # ---------- generate strained structure
    while True:
        # ------ strain matrix
        strain_matrix = np.empty([3, 3])
        for i in range(3):
            for j in range(3):
                if i <= j:
                    if i == j:
                        strain_matrix[i][j] = 1.0 + np.random.normal(loc=0.0, scale=sigma_st)
                    else:
                        strain_matrix[i][j] = np.random.normal(loc=0.0, scale=sigma_st)/2.0
                        strain_matrix[j][i] = strain_matrix[i][j]
        # ------ strained lattice
        strained_lattice = np.dot(strain_matrix, lat_mat).T
        # ------ child
        child = Structure(strained_lattice, child.species, child.frac_coords)
        # ------ scale lattice
        child.scale_lattice(parent_A.volume)
        # ------ check distance
        success, mindist_ij, dist = check_distance(child, atype, mindist)
        if success:
            child = sort_by_atype(child, atype)
            return child
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in strain: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                return None    # change parent

