import configparser
from .base import BaseReader


class EAReader(BaseReader):
    """
    Reader for [EA] section
    """
    def read(self):
        # ---------- n_pop
        self.rin.n_pop = self.config.getint('EA', 'n_pop')
        if self.rin.n_pop < 1:
            raise ValueError('n_pop must be positive int')

        # ---------- n_crsov
        self.rin.n_crsov = self.config.getint('EA', 'n_crsov')
        if self.rin.n_crsov < 0:
            raise ValueError('n_crsov must be non-negative int')

        # ---------- n_perm
        self.rin.n_perm = self.config.getint('EA', 'n_perm')
        if self.rin.n_perm < 0:
            raise ValueError('n_perm must be non-negative int')
        if self.rin.n_perm != 0 and len(self.rin.atype) == 1:
            raise ValueError('When the number of atom type is 1, n_perm must be 0')

        # ---------- n_strain
        self.rin.n_strain = self.config.getint('EA', 'n_strain')
        if self.rin.n_strain < 0:
            raise ValueError('n_strain must be non-negative int')

        # ---------- n_rand
        self.rin.n_rand = self.config.getint('EA', 'n_rand')
        if self.rin.n_rand < 0:
            raise ValueError('n_rand must be non-negative int')

        # ---------- check n_pop for EA
        if self.rin.algo == 'EA' and self.rin.struc_mode not in ['mol', 'mol_bs']:
            if self.rin.n_crsov + self.rin.n_perm + self.rin.n_strain + self.rin.n_rand != self.rin.n_pop:
                raise ValueError('n_crsov + n_perm + n_strain + n_rand must be n_pop')

        # ---------- n_elite
        self.rin.n_elite = self.config.getint('EA', 'n_elite')
        if self.rin.n_elite < 0:
            raise ValueError('n_elite must be non-negative int')

        # ---------- fit_reverse
        try:
            self.rin.fit_reverse = self.config.getboolean('EA', 'fit_reverse')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.fit_reverse = False

        # ---------- n_fittest
        try:
            self.rin.n_fittest = self.config.getint('EA', 'n_fittest')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.n_fittest = 0
        if self.rin.n_fittest < 0:
            raise ValueError('n_fittest must be non-negative int')

        # ---------- slct_func
        self.rin.slct_func = self.config.get('EA', 'slct_func')
        if self.rin.slct_func not in ['TNM', 'RLT']:
            raise ValueError('slct_func must be TNM or RLT')
        # ------ TNM
        if self.rin.slct_func == 'TNM':
            # -- t_size
            try:
                self.rin.t_size = self.config.getint('EA', 't_size')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.t_size = 3
            if self.rin.t_size < 2:
                raise ValueError('t_size must be greater than or equal to 2')
        # ------ RLT
        elif self.rin.slct_func == 'RLT':
            # -- a_rlt
            try:
                self.rin.a_rlt = self.config.getfloat('EA', 'a_rlt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.a_rlt = 10.0
            # -- b_rlt
            try:
                self.rin.b_rlt = self.config.getfloat('EA', 'b_rlt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.b_rlt = 1.0
            # -- check
            if not 0 < self.rin.b_rlt < self.rin.a_rlt:
                raise ValueError('must be 0 < b_rlt < a_rlt')

        # ---------- crossover
        if self.rin.n_crsov > 0:
            # ------ crs_lat
            try:
                self.rin.crs_lat = self.config.get('EA', 'crs_lat')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.crs_lat = 'random'
            if self.rin.crs_lat not in ['equal', 'random']:
                raise ValueError('crs_lat must be equal or random')
            # ------ nat_diff_tole
            try:
                self.rin.nat_diff_tole = self.config.getint('EA', 'nat_diff_tole')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.nat_diff_tole = 4
            if self.rin.nat_diff_tole < 0:
                raise ValueError('nat_diff_tole must be non-negative int')

        # ---------- permutation
        if self.rin.n_perm > 0:
            # ------ ntimes
            try:
                self.rin.ntimes = self.config.getint('EA', 'ntimes')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.ntimes = 1
            if self.rin.ntimes < 1:
                raise ValueError('ntimes must be positive int')

        # ---------- strain
        if self.rin.n_strain > 0:
            # ------ sigma_st
            try:
                self.rin.sigma_st = self.config.getfloat('EA', 'sigma_st')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.sigma_st = 0.5
            if self.rin.sigma_st <= 0:
                raise ValueError('sigma_st must be positive float')

        # ---------- maxcnt_ea
        try:
            self.rin.maxcnt_ea = self.config.getint('EA', 'maxcnt_ea')
        except configparser.NoOptionError:
            self.rin.maxcnt_ea = 50

        # ---------- maxgen_ea
        try:
            self.rin.maxgen_ea = self.config.getint('EA', 'maxgen_ea')
        except configparser.NoOptionError:
            self.rin.maxgen_ea = 0
        if self.rin.maxgen_ea < 0:
            raise ValueError('maxgen_ea must be non-negative int')

        # ---------- emin_ea, emax_ea
        try:
            self.rin.emin_ea = self.config.getfloat('EA', 'emin_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.emin_ea = None
        try:
            self.rin.emax_ea = self.config.getfloat('EA', 'emax_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.emax_ea = None
        if self.rin.emax_ea is not None and self.rin.emin_ea is not None:
            if self.rin.emin_ea > self.rin.emax_ea:
                raise ValueError('emax_ea < emin_ea, check emax_ea and emin_ea')

        # ---------- EA-vc
        if self.rin.algo == 'EA-vc':
            self._read_ea_vc()

        # ---------- mol or mol_bs
        if self.rin.struc_mode in ['mol', 'mol_bs']:
            self._read_ea_mol()


    def _read_ea_vc(self):
        # ---------- n_add
        self.rin.n_add = self.config.getint('EA', 'n_add')
        if self.rin.n_add < 0:
            raise ValueError('n_add must be non-negative int')

        # ---------- add_max
        try:
            self.rin.add_max = self.config.getint('EA', 'add_max')
            if self.rin.add_max <= 0:
                raise ValueError('add_max must be non-negative int')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.add_max = 3

        # ---------- n_elim
        self.rin.n_elim = self.config.getint('EA', 'n_elim')
        if self.rin.n_elim < 0:
            raise ValueError('n_elim must be non-negative int')

        # ---------- elim_max
        try:
            self.rin.elim_max = self.config.getint('EA', 'elim_max')
            if self.rin.elim_max <= 0:
                raise ValueError('elim_max must be non-negative int')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.elim_max = 3

        # ---------- n_subs
        self.rin.n_subs = self.config.getint('EA', 'n_subs')
        if self.rin.n_subs < 0:
            raise ValueError('n_subs must be non-negative int')
        if self.rin.n_subs != 0:
            diff_nat = [u - l for u, l in zip(self.rin.ul_nat, self.rin.ll_nat)]
            non_zero_count = len([x for x in diff_nat if x != 0])
            if non_zero_count < 2:
                raise ValueError('n_subs: at least two atomic types are variable.')

        # ---------- subs_max
        try:
            self.rin.subs_max = self.config.getint('EA', 'subs_max')
            if self.rin.subs_max <= 0:
                raise ValueError('subs_max must be non-negative int')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.subs_max = 3

        # ---------- target
        self.rin.target = self.config.get('EA', 'target')
        if self.rin.target not in ['random']:
            raise ValueError('target must be random for now')

        # ---------- check n_pop for EA-vc
        if self.rin.n_crsov + self.rin.n_perm + self.rin.n_strain + self.rin.n_rand \
            + self.rin.n_add + self.rin.n_elim + self.rin.n_subs != self.rin.n_pop:
            raise ValueError('n_crsov + n_perm + n_strain + n_rand'
                                ' + n_add + n_elim + n_subs must be n_pop')

        # ---------- end_point
        self.rin.end_point = self.config.get('EA', 'end_point')
        self.rin.end_point = tuple([float(x) for x in self.rin.end_point.split()])
        if not len(self.rin.end_point) == len(self.rin.atype):
            raise ValueError('len(end_point) == len(atype), check end_point')

        # ---------- cgen
        try:
            self.rin.cgen = self.config.getint('EA', 'cgen')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.cgen = None

        # ---------- show_max
        try:
            self.rin.show_max = self.config.getfloat('EA', 'show_max')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.show_max = 0.2

        # ---------- label_stable
        try:
            self.rin.label_stable = self.config.getboolean('EA', 'label_stable')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.label_stable = True

        # ---------- vmax
        try:
            self.rin.vmax = self.config.getfloat('EA', 'vmax')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.vmax = 0.2

        # ---------- bottom_margin
        try:
            self.rin.bottom_margin = self.config.getfloat('EA', 'bottom_margin')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.bottom_margin = 0.02

        # ---------- fig_format
        try:
            self.rin.fig_format = self.config.get('EA', 'fig_format')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.fig_format = 'svg'
        if self.rin.fig_format not in ['svg', 'png', 'pdf']:
            raise ValueError('fig_format must be svg, png, or pdf')


    def _read_ea_mol(self):
        # ---------- n_rotation
        self.rin.n_rotation = self.config.getint('EA', 'n_rotation')
        if self.rin.n_rotation < 0:
            raise ValueError('n_rotation must be non-negative int')

        # ---------- check n_pop for mol or mol_bs
        if self.rin.n_crsov + self.rin.n_perm + self.rin.n_strain + self.rin.n_rand \
            + self.rin.n_rotation != self.rin.n_pop:
            raise ValueError('n_crsov + n_perm + n_strain + n_rand + n_rotation must be n_pop')

        # ---------- mindist_mol_ea
        self.rin.mindist_mol_ea = []
        for i in range(len(self.rin.mol_file)):
            tmp = self.config.get('EA', f'mindist_mol_ea_{i+1}')
            tmp = tuple([float(x) for x in tmp.split()])
            if not len(tmp) == len(self.rin.mol_file):
                raise ValueError(f'not len(mindist_mol_ea_{i+1}) == len(mol_file)')
            self.rin.mindist_mol_ea.append(tmp)
        self.rin.mindist_mol_ea = tuple(self.rin.mindist_mol_ea)    # list to tuple
        # ------ check symmetric matrix
        for i in range(len(self.rin.mindist_mol_ea)):
            for j in range(len(self.rin.mindist_mol_ea)):
                if i < j:
                    if not self.rin.mindist_mol_ea[i][j] == self.rin.mindist_mol_ea[j][i]:
                        raise ValueError(f'mindist_mol_ea is not symmetric. ({i}, {j}) -->'
                                            f' {self.rin.mindist_mol_ea[i][j]}, ({j}, {i}) --> {self.rin.mindist_mol_ea[j][i]}')

        # ---------- rot_max_angle
        try:
            self.rin.rot_max_angle = self.config.getint('EA', 'rot_max_angle')
        except configparser.NoOptionError:
            self.rin.rot_max_angle = 360
        if self.rin.rot_max_angle <= 0:
            raise ValueError('rot_max_angle must be positive int')

        # ---------- protect_mol_struc
        try:
            self.rin.protect_mol_struc = self.config.getboolean('EA', 'protect_mol_struc')
        except configparser.NoOptionError:
            self.rin.protect_mol_struc = False