'''
Structure generation by evolutionary algorithm
'''

from logging import getLogger

from ...IO import read_input as rin
from ...util.struc_util import out_poscar


logger = getLogger('cryspy')

class EA_generation:
    '''
    generate structures by evolutionary algorithm

    # ---------- args
    sp: instance of Select_parents class
    id_start (int or None): starting id
                            if None, id_start = max(self.fitness.keys()) + 1
    init_pos_path (str or None): if not None, structure data in POSCAR format
                                     is appended to init_pos_path

    # ---------- instance methods
    self.gen_crossover(self, co)

    self.gen_permutation(self, pm)

    self.gen_strain(self, st)
    '''

    def __init__(self, sp, id_start=None, init_pos_path=None):
        # ---------- check args
        self.sp = sp
        # ------ id_offset
        if id_start is None:
            self.cid = max(sp.fitness.keys()) + 1
        elif isinstance(id_start, int):
            if id_start < (max(sp.fitness.keys()) + 1):
                logger.error('id_start is already included'
                                 ' structure ID of the data')
            else:
                self.cid = id_start
        else:
            logger.error('id_start must be int or None')
        # ------ init_pos_path
        if init_pos_path is not None:
            if isinstance(init_pos_path, str):
                self.init_pos_path = init_pos_path
            else:
                logger.error('init_pos_path must be str or None')
        # ---------- initialize data
        self.offspring = {}    # structure data
        self.offspring_mol_id = {}
        self.parents = {}    # tuple of parents ID
        self.operation = {}

    def gen_crossover(self, co, struc_mol_id, molecular=False):
        '''
        generate structures by crossover

        # ---------- args
        co: instance of Crossover class
        '''
        # ---------- generate structures by crossover
        struc_cnt = 0
        while struc_cnt < rin.n_crsov:
            # ------ select parents
            pid_A, pid_B = self.sp.get_parents(n_parent=2)
            # ------ generate child
            if molecular:
                child, mol_id = co.gen_child_mol(self.sp.struc_data[pid_A], self.sp.struc_data[pid_B],
                                                 struc_mol_id[pid_A], struc_mol_id[pid_B])
            else:
                child = co.gen_child(self.sp.struc_data[pid_A],
                                     self.sp.struc_data[pid_B])
            # ------ success
            if child is not None:
                self.offspring[self.cid] = child
                if molecular:
                    self.offspring_mol_id[self.cid] = mol_id
                self.parents[self.cid] = (pid_A, pid_B)
                self.operation[self.cid] = 'crossover'
                try:
                    spg_sym, spg_num = child.get_space_group_info(
                        symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                logger.info(f'Structure ID {self.cid:>6} was generated'
                      f' from {pid_A:>6} and {pid_B:>6} by crossover.'
                      f' Space group: {spg_num:>3} {spg_sym}')
                if self.init_pos_path is not None:
                    out_poscar(child, self.cid, self.init_pos_path)
                self.cid += 1
                struc_cnt += 1

    def gen_permutation(self, pm, struc_mol_id, molecular=False):
        '''
        generate structures by permutation

        # ---------- args
        pm: instance of Permutation class
        '''
        # ---------- generate structures by permutation
        struc_cnt = 0
        while struc_cnt < rin.n_perm:
            # ------ select parents
            pid, = self.sp.get_parents(n_parent=1)    # comma for list[0]
            # ------ generate child
            if molecular:
                child, mol_id = pm.gen_child_mol(self.sp.struc_data[pid], struc_mol_id[pid])
            else:
                child = pm.gen_child(self.sp.struc_data[pid])
            # ------ success
            if child is not None:
                self.offspring[self.cid] = child
                if molecular:
                    self.offspring_mol_id[self.cid] = mol_id
                self.parents[self.cid] = (pid, )    # tuple
                self.operation[self.cid] = 'permutation'
                try:
                    spg_sym, spg_num = child.get_space_group_info(
                        symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                logger.info(f'Structure ID {self.cid:>6} was generated'
                      f' from {pid:>6} by permutation.'
                      f' Space group: {spg_num:>3} {spg_sym}')
                if self.init_pos_path is not None:
                    out_poscar(child, self.cid, self.init_pos_path)
                self.cid += 1
                struc_cnt += 1

    def gen_strain(self, st, struc_mol_id, molecular=False, protect_mol_struc=True):
        '''
        generate structures by strain

        # ---------- args
        st: instance of Strain class
        '''
        # ---------- generate structures by strain
        struc_cnt = 0
        while struc_cnt < rin.n_strain:
            # ------ select parents
            pid, = self.sp.get_parents(n_parent=1)    # comma for list[0]
            # ------ generate child
            if molecular:
                if protect_mol_struc:
                    child, mol_id = st.gen_child_mol(self.sp.struc_data[pid], struc_mol_id[pid])
                else:
                    child = st.gen_child(self.sp.struc_data[pid])
                    mol_id = struc_mol_id[pid]
            else:
                child = st.gen_child(self.sp.struc_data[pid])
            # ------ success
            if child is not None:
                self.offspring[self.cid] = child
                if molecular:
                    self.offspring_mol_id[self.cid] = mol_id
                self.parents[self.cid] = (pid, )    # tuple
                self.operation[self.cid] = 'strain'
                try:
                    spg_sym, spg_num = child.get_space_group_info(
                        symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                logger.info(f'Structure ID {self.cid:>6} was generated'
                      f' from {pid:>6} by strain.'
                      f' Space group: {spg_num:>3} {spg_sym}')
                if self.init_pos_path is not None:
                    out_poscar(child, self.cid, self.init_pos_path)
                self.cid += 1
                struc_cnt += 1

    def gen_rotation(self, struc_mol_id, rot):
        '''
        generate structures by rotation
        only for mol

        # ---------- args
        st: instance of rotation class
        '''
        # ---------- generate structures by rotation
        struc_cnt = 0
        while struc_cnt < rin.n_rotation:
            # ------ select parents
            pid, = self.sp.get_parents(n_parent=1)    # comma for list[0]
            # ------ generate child
            child, mol_id = rot.gen_child(self.sp.struc_data[pid], struc_mol_id[pid])
            # ------ success
            if child is not None:
                self.offspring[self.cid] = child
                self.offspring_mol_id[self.cid] = mol_id
                self.parents[self.cid] = (pid, )    # tuple
                self.operation[self.cid] = 'rotation'
                try:
                    spg_sym, spg_num = child.get_space_group_info(symprec=rin.symprec)
                except TypeError:
                    spg_num = 0
                    spg_sym = None
                logger.info(f'Structure ID {self.cid:>6} was generated'
                      f' from {pid:>6} by rotation.'
                      f' Space group: {spg_num:>3} {spg_sym}')
                if self.init_pos_path is not None:
                    out_poscar(child, self.cid, self.init_pos_path)
                self.cid += 1
                struc_cnt += 1
