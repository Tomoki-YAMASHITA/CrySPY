import configparser
from dataclasses import dataclass, field
from logging import getLogger
import os


logger = getLogger('cryspy')


@dataclass
class ReadInput:
    '''
    Read cryspy.in
        Here, ignore type hint for None
    '''
    # ---------- basic section
    algo: str = field(default=None)
    calc_code: str = field(default=None)
    tot_struc: int = field(default=None)
    nstage: int = field(default=None)
    njob: int = field(default=None)
    jobcmd: str = field(default=None)
    jobfile: str = field(default=None)

    # ---------- structure section
    struc_mode: str = field(default=None)
    atype: tuple = field(default=None)
    nat: tuple = field(default=None)
    mindist: tuple = field(default=None)
    mindist_factor: float = field(default=None)
    vol_factor: float = field(default=None)
    vol_mu: float = field(default=None)
    vol_sigma: float = field(default=None)
    symprec: float = field(default=None)
    spgnum: str = field(default=None)    # ignore type hint here. 'all', tuple or 0
    use_find_wy: bool = field(default=None)
    # ------ EA-vc
    ll_nat: tuple = field(default=None)
    ul_nat: tuple = field(default=None)
    charge: tuple = field(default=None)
    cn_nmax: int = field(default=None)
    # ------ mol or mol_bs
    mol_file: tuple = field(default=None)
    nmol: tuple = field(default=None)
    timeout_mol: float = field(default=None)
    rot_mol: str = field(default=None)
    nrot: int = field(default=None)
    mindist_mol_bs: tuple = field(default=None)
    mindist_mol_bs_factor: float = field(default=None)
    # ------ find_wy or spgnum = 0
    fwpath: str = field(default=None)
    minlen: float = field(default=None)
    maxlen: float = field(default=None)
    dangle: float = field(default=None)
    maxcnt: int = field(default=None)

    # ---------- option section
    stop_chkpt: int = field(default=None)
    load_struc_flag: bool = field(default=None)
    stop_next_struc: bool = field(default=None)
    recalc: tuple = field(default=None)
    append_struc_ea: bool = field(default=None)
    energy_step_flag: bool = field(default=None)
    struc_step_flag: bool = field(default=None)
    force_step_flag: bool = field(default=None)
    stress_step_flag: bool = field(default=None)

    # ---------- BO section
    nselect_bo: int = field(default=None)
    score: str = field(default=None)
    num_rand_basis: int = field(default=None)
    cdev: float = field(default=None)
    dscrpt: str = field(default=None)
    fp_rmin: float = field(default=None)
    fp_rmax: float = field(default=None)
    fp_npoints: int = field(default=None)
    fp_sigma: float = field(default=None)
    max_select_bo: int = field(default=None)
    manual_select_bo: tuple = field(default=None)
    emin_bo: float = field(default=None)
    emax_bo: float = field(default=None)

    # ---------- LAQA section
    nselect_laqa: int = field(default=None)
    wf: float = field(default=None)
    ws: float = field(default=None)

    # ---------- EA section
    n_pop: int = field(default=None)
    n_crsov: int = field(default=None)
    n_perm: int = field(default=None)
    n_strain: int = field(default=None)
    n_rand: int = field(default=None)
    n_elite: int = field(default=None)
    fit_reverse: bool = field(default=None)
    n_fittest: int = field(default=None)
    slct_func: str = field(default=None)
    t_size: int = field(default=None)
    a_rlt: float = field(default=None)
    b_rlt: float = field(default=None)
    crs_lat: str = field(default=None)
    nat_diff_tole: int = field(default=None)
    ntimes: int = field(default=None)
    sigma_st: float = field(default=None)
    maxcnt_ea: int = field(default=None)
    maxgen_ea: int = field(default=None)
    emin_ea: float = field(default=None)
    emax_ea: float = field(default=None)
    n_add: int = field(default=None)
    n_elim: int = field(default=None)
    n_subs: int = field(default=None)
    target: str = field(default=None)
    end_point: tuple = field(default=None)
    show_max: float = field(default=None)
    lable_stable: bool = field(default=None)
    vmax: float = field(default=None)
    bottom_margin: float = field(default=None)
    n_rotation: int = field(default=None)          # not implemented yet, for EA mol
    mindist_mol_ea: tuple = field(default=None)    # not implemented yet, for EA mol
    rot_max_angle: float = field(default=None)     # not implemented yet, for EA mol
    protect_mol_struc: bool = field(default=None)  # not implemented yet, for EA mol

    # ---------- common in VASP, QE, OMX
    kpt_flag: bool = field(default=None)
    kppvol: int = field(default=None)
    force_gamma: bool = field(default=None)

    # ---------- VASP section
    # no inputs other than common

    # ---------- QE section
    qe_infile: str = field(default=None)
    qe_outfile: str = field(default=None)
    pv_term: bool = field(default=None)

    # ---------- OMX section
    OMX_infile: str = field(default=None)
    OMX_outfile: str = field(default=None)
    upSpin: dict = field(default=None)
    downSpin: dict = field(default=None)

    # ---------- soiap section
    soiap_infile: str = field(default=None)
    soiap_outfile: str = field(default=None)
    soiap_cif: str = field(default=None)

    # ---------- lammps section
    lammps_infile: str = field(default=None)
    lammps_outfile: str = field(default=None)
    lammps_potential: str = field(default=None)
    lammps_data: str = field(default=None)

    # ---------- ASE section
    ase_python: str = field(default=None)

    def __post_init__(self):
        self._read_cryspyin()    # self.config is created
        self._read_basic()       # read basic section in cryspy.in
        self._read_structure()   # read structure section in cryspy.in
        self._read_option()      # read option section in cryspy.in, do this prior to EA for append_struc_ea
        if self.algo == 'BO':
            self._read_bo()      # read BO section in cryspy.in
        if self.algo == 'LAQA':
            self._read_laqa()    # read LAQA section in cryspy.in
        if self.algo in ['EA', 'EA-vc'] or self.append_struc_ea:
            self._read_ea()      # read EA section in cryspy.in
        if self.calc_code == 'VASP':
            self._read_vasp()    # read VASP section in cryspy.in
        if self.calc_code == 'QE':
            self._read_qe()      # read QE section in cryspy.in
        if self.calc_code == 'OMX':
            self._read_omx()     # read OMX section in cryspy.in
        if self.calc_code == 'soiap':
            self._read_soiap()   # read soiap section in cryspy.in
        if self.calc_code == 'LAMMPS':
            self._read_lammps()  # read LAMMPS section in cryspy.in
        if self.calc_code == 'ASE':
            self._read_ase()     # read ASE section in cryspy.in

    def _read_cryspyin(self):
        filename = 'cryspy.in'
        if not os.path.isfile(filename):
            raise FileNotFoundError(f'{filename} not found in {os.getcwd()}')
        self.config = configparser.ConfigParser()
        self.config.read(filename)

    def _read_basic(self):
        # ---------- algo
        self.algo = self.config.get('basic', 'algo')
        if self.algo not in ['RS', 'BO', 'LAQA', 'EA', 'EA-vc']:
            raise ValueError('algo must be RS, BO, LAQA, EA or EA-vc')
        # ---------- calc_code
        self.calc_code = self.config.get('basic', 'calc_code')
        if self.calc_code not in ['VASP', 'QE', 'soiap', 'LAMMPS', 'OMX', 'ASE']:
            raise ValueError('calc_code must be VASP, QE, OMX, soiap, LAMMPS, or ASE')
        if self.algo == 'LAQA':
            if self.calc_code not in ['VASP', 'QE', 'soiap']:
                raise ValueError('LAQA: only VASP, QE, or soiap for now')
        # ---------- tot_struc
        if self.algo not in ['EA', 'EA-vc']:
            self.tot_struc = self.config.getint('basic', 'tot_struc')
            if self.tot_struc < 1:
                raise ValueError('tot_struc < 1, check tot_struc')
        # ---------- nstage
        self.nstage = self.config.getint('basic', 'nstage')
        if self.nstage < 1:
            raise ValueError('nstage < 1, check nstage')
        if self.algo == 'LAQA':
            if not self.nstage == 1:
                raise ValueError('nstage must be 1 in LAQA')
        # ---------- njob
        self.njob = self.config.getint('basic', 'njob')
        if self.njob < 1:
            raise ValueError('njob < 1, check njob')
        # ---------- jobcmd
        self.jobcmd = self.config.get('basic', 'jobcmd')
        # ---------- jobfile
        self.jobfile = self.config.get('basic', 'jobfile')

    def _read_structure(self):
        # ---------- struc_mode
        try:
            self.struc_mode = self.config.get('structure', 'struc_mode')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.struc_mode = 'crystal'
        if self.struc_mode not in ['crystal', 'mol', 'mol_bs']:
            raise ValueError('struc_mode must be crystal, mol, or  mol_bs')
        #if self.struc_mode not in ['crystal', 'mol', 'mol_bs', 'host']:
        #    raise ValueError('struc_mode must be crystal, mol, mol_bs, or host')
        if self.algo in ['EA', 'EA-vc'] and self.struc_mode in ['mol', 'mol_bs']:
                raise ValueError('currently, molecular structure is not allowed in EA or EA-vc')
        # ---------- atype
        self.atype = self.config.get('structure', 'atype')
        self.atype = tuple([a for a in self.atype.split()])    # str --> list --> tuple
        if self.algo == 'EA-vc':
            if len(self.atype) < 2:
                raise ValueError('EA-vc: atype must have at least 2 elements')
        # ---------- nat
        if not self.algo == 'EA-vc':
            self.nat = self.config.get('structure', 'nat')
            self.nat = tuple([int(x) for x in self.nat.split()])    # str --> int --> list --> tuple
            if not len(self.nat) == len(self.atype):
                raise ValueError('not len(nat) == len(atype), check atype and nat')
        # ---------- mindist
        try:
            self.mindist = []
            for i in range(len(self.atype)):
                tmp = self.config.get('structure', f'mindist_{i+1}')
                tmp = tuple([float(x) for x in tmp.split()])    # str --> float --> list --> tuple
                if not len(tmp) == len(self.atype):
                    raise ValueError(f'not len(mindist_{i+1}) == len(atype)')
                self.mindist.append(tmp)
            self.mindist = tuple(self.mindist)    # list --> tuple
            # ------ check symmetric matrix
            for i in range(len(self.mindist)):
                for j in range(len(self.mindist)):
                    if i < j:
                        if not self.mindist[i][j] == self.mindist[j][i]:
                            raise ValueError(f'mindist is not symmetric. ({i}, {j}) -->'
                                            f' {self.mindist[i][j]}, ({j}, {i}) --> {self.mindist[j][i]}')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.mindist = None
        # ---------- mindist_factor
        try:
            self.mindist_factor = self.config.getfloat('structure', 'mindist_factor')
            if self.mindist_factor <= 0.0:
                raise ValueError('mindist_factor must be positive')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.mindist_factor = 1.0
        # ---------- vol_factor
        try:
            self.vol_factor = self.config.getfloat('structure', 'vol_factor')
            if self.vol_factor <= 0.0:
                raise ValueError('vol_factor must be positive')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.vol_factor = 1.1
        # ---------- vol_mu
        try:
            self.vol_mu = self.config.getfloat('structure', 'vol_mu')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.vol_mu = None
        if self.vol_mu is not None:
            if self.vol_mu <= 0.0:
                raise ValueError('vol_mu must be positive float')
        # ---------- vol_sigma
        try:
            self.vol_sigma = self.config.getfloat('structure', 'vol_sigma')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.vol_sigma = None
        if self.vol_mu is not None and self.vol_sigma is None:
            raise ValueError('vol_mu is set, but vol_sigma is not set')
        if self.vol_sigma is not None:
            if self.vol_sigma < 0.0:
                raise ValueError('vol_sigma must not be negative')
        # ---------- symprec
        try:
            self.symprec = self.config.getfloat('structure', 'symprec')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.symprec = 0.01
        if self.symprec < 0.0:
            raise ValueError('symprec must be positive float')
        # ---------- spgnum
        try:
            self.spgnum = self.config.get('structure', 'spgnum')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.spgnum = 'all'
        if self.spgnum == '0':    # str here
            if self.struc_mode in ['mol', 'mol_bs']:
                raise ValueError('spgnum = 0 is not allow when struc_mode is mol or mol_bs')
            self.spgnum = 0    # int here
        elif self.spgnum == 'all':
            pass
        else:
            self._spgtuple()    # self.spgnum: e.g. 3 15-18 100 --> (3, 15, 16, 17, 18, 100)
        # ---------- use_find_wy
        try:
            self.use_find_wy = self.config.getboolean('structure', 'use_find_wy')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.use_find_wy = False
        if self.use_find_wy:
            if not self.struc_mode == 'crystal':
                raise ValueError('find_wy can be use if struc_mode is crystal')
        # ---------- EA-vc
        if self.algo =='EA-vc':
            # ------ ll_nat, ul_nat
            self.ll_nat = self.config.get('structure', 'll_nat')
            self.ll_nat = tuple([int(x) for x in self.ll_nat.split()])    # str --> int --> list --> tuple
            self.ul_nat = self.config.get('structure', 'ul_nat')
            self.ul_nat = tuple([int(x) for x in self.ul_nat.split()])    # str --> int --> list --> tuple
            if not len(self.atype) == len(self.ll_nat) == len(self.ul_nat):
                raise ValueError('not len(atype) == len(ll_nat) == len(ul_nat), check ll_nat and ul_nat')
            for i in range(len(self.ll_nat)):
                if not 0 <= self.ll_nat[i] <= self.ul_nat[i]:
                    raise ValueError(f'not 1 <= ll_nat[{i}] <= ul_nat[{i}], check ll_nat and ul_nat')
            # ------ charge
            try:
                self.charge = self.config.get('structure', 'charge')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.charge = None
            if self.charge is not None:
                self.charge = tuple([int(x) for x in self.charge.split()])    # str --> tuple
                if not len(self.charge) == len(self.atype):
                    raise ValueError('not len(charge) == len(atype), check charge')
                if any(x > 0 for x in self.charge) and any(x < 0 for x in self.charge):
                    pass
                else:
                    raise ValueError('charge must have both positive and negative integers')
            # ------ cn_nmax
            try:
                self.cn_nmax = self.config.getint('structure', 'cn_nmax')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.cn_nmax = None
            if self.cn_nmax is not None:
                if self.cn_nmax < 2:
                    raise ValueError('cn_nmax must be more than 1')
            if self.charge is not None and self.cn_nmax is None:
                self.cn_nmax = 3
        # ---------- mol or mol_bs
        if self.struc_mode in ['mol', 'mol_bs']:
            # ------ mol_file
            self.mol_file = self.config.get('structure', 'mol_file')
            self.mol_file = tuple([a for a in self.mol_file.split()])    # str --> list --> tuple
            # ------ nmol
            self.nmol = self.config.get('structure', 'nmol')
            self.nmol = tuple([int(x) for x in self.nmol.split()])    # str --> int --> list --> tuple
            if not len(self.mol_file) == len(self.nmol):
                raise ValueError('not len(mol_file) == len(nmol)')
            # ------ timeout_mol
            try:
                self.timeout_mol = self.config.getfloat('structure', 'timeout_mol')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.timeout_mol = None
            if self.timeout_mol is not None:
                if self.timeout_mol <= 0:
                    raise ValueError('timeout_mol must be None or positive')
        # ---------- mol_bs
        if self.struc_mode == 'mol_bs':
            # ------ rot_mol
            try:
                self.rot_mol = self.config.get('structure', 'rot_mol')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rot_mol = 'random_wyckoff'
            if self.rot_mol not in ['random', 'random_mol', 'random_wyckoff']:
                raise ValueError('rot_mol must be random, random_mol, or random_wyckoff')
            # ------ nrot
            try:
                self.nrot = self.config.getint('structure', 'nrot')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.nrot = 20
            if self.nrot < 1:
                raise ValueError('nrot < 1, check nrot')
            # ------ mindist_mol_bs
            try:
                self.mindist_mol_bs = []
                for i in range(len(self.mol_file)):
                    tmp = self.config.get('structure', f'mindist_mol_bs_{i+1}')
                    tmp = tuple([float(x) for x in tmp.split()])    # str --> float --> list --> tuple
                    if not len(tmp) == len(self.mol_file):
                        raise ValueError(f'not len(mindist_mol_bs_{i+1}) == len(mol_file)')
                    self.mindist_mol_bs.append(tmp)
                self.mindist_mol_bs = tuple(self.mindist_mol_bs)    # list --> tuple
                # -- check symmetric matrix
                for i in range(len(self.mindist_mol_bs)):
                    for j in range(len(self.mindist_mol_bs)):
                        if i < j:
                            if not self.mindist_mol_bs[i][j] == self.mindist_mol_bs[j][i]:
                                raise ValueError(f'mindist_mol_bs is not symmetric. ({i}, {j}) -->'
                                                f' {self.mindist_mol_bs[i][j]}, ({j}, {i}) --> {self.mindist_mol_bs[j][i]}')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.mindist_mol_bs = None
            # ------ mindist_mol_bs_factor
            try:
                self.mindist_mol_bs_factor = self.config.getfloat('structure', 'mindist_mol_bs_factor')
                if self.mindist_mol_bs_factor <= 0.0:
                    raise ValueError('mindist_mol_bs_factor must be positive')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.mindist_mol_bs_factor = 1.0
        # ---------- find_wy or spgnum = 0
        if self.spgnum == 0 or self.use_find_wy:
            # ------ fwpath
            try:
                self.fwpath = self.config.get('structure', 'fwpath')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.fwpath = None
            # ------ minlen
            self.minlen = self.config.getfloat('structure', 'minlen')
            if self.minlen <= 0.0:
                raise ValueError('minlen must be positive')
            # ------ maxlen
            self.maxlen = self.config.getfloat('structure', 'maxlen')
            if self.minlen > self.maxlen:
                raise ValueError('minlen > maxlen')
            # ------ dangle
            self.dangle = self.config.getfloat('structure', 'dangle')
            if self.dangle <= 0.0:
                raise ValueError('dangle must be positive')
            # ------ maxcnt
            try:
                self.maxcnt = self.config.getint('structure', 'maxcnt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.maxcnt = 50
            if self.maxcnt < 1:
                raise ValueError('maxcnt must be positive int')

    def _spgtuple(self):
        tmpspg = []
        for c in self.spgnum.split():
            if '-' in c:
                if not len(c.split('-')) == 2:
                    raise ValueError('Wrong input in spgnum.')
                istart = int(c.split('-')[0])
                iend = int(c.split('-')[1])+1
                if istart < 0 or 230 < istart:
                    raise ValueError('spgnum must be 1 -- 230')
                if iend < 0 or 231 < iend:
                    raise ValueError('spgnum must be 1 -- 230')
                for i in range(istart, iend):
                    if i not in tmpspg:
                        tmpspg.append(i)
            else:
                if int(c) < 0 or 230 < int(c):
                    raise ValueError('spgnum must be 1 -- 230')
                if not int(c) in tmpspg:
                    tmpspg += [int(c)]
        self.spgnum = tuple(tmpspg)

    def _read_option(self):
        # ---------- stop_chkpt
        try:
            self.stop_chkpt = self.config.getint('option', 'stop_chkpt')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.stop_chkpt = 0
        # ---------- load_struc_flag
        try:
            self.load_struc_flag = self.config.getboolean('option', 'load_struc_flag')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.load_struc_flag = False
        # ---------- stop_next_struc
        try:
            self.stop_next_struc = self.config.getboolean('option', 'stop_next_struc')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.stop_next_struc = False
        # ---------- recalc
        try:
            self.recalc = self.config.get('option', 'recalc')
            self.recalc = tuple([int(x) for x in self.recalc.split()])    # str --> int --> list --> tuple
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.recalc = None
        if self.recalc is not None:
            if self.algo not in ['EA', 'EA-vc']:
                for cid in self.recalc:
                    if not 0 <= cid < self.tot_struc:
                        raise ValueError('recalc must be non-negative int and less than tot_struc')
        # ---------- append_struc_ea
        try:
            self.append_struc_ea = self.config.getboolean('option', 'append_struc_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.append_struc_ea = False
        # ---------- energy_step_flag
        try:
            self.energy_step_flag = self.config.getboolean('option', 'energy_step_flag')
            if self.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('energy_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.energy_step_flag = False
        try:
            self.struc_step_flag = self.config.getboolean('option', 'struc_step_flag')
            if self.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('struc_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.struc_step_flag = False
        try:
            self.force_step_flag = self.config.getboolean('option', 'force_step_flag')
            if self.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('force_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.force_step_flag = False
        if self.algo == 'LAQA':
            self.force_step_flag = True
        try:
            self.stress_step_flag = self.config.getboolean('option', 'stress_step_flag')
            if self.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('stress_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.stress_step_flag = False
        if self.algo == 'LAQA':
            self.stress_step_flag = True

    def _read_bo(self):
        # ---------- check physbo
        try:
            import physbo
        except ModuleNotFoundError:
            raise ModuleNotFoundError('PHYSBO is required for BO. --> pip3 install physbo')
        # ---------- nselect_bo
        self.nselect_bo = self.config.getint('BO', 'nselect_bo')
        if self.nselect_bo < 1:
            raise ValueError('nselect_bo < 1, check nselect_bo')
        elif self.tot_struc < self.nselect_bo:
            raise ValueError('tot_struc < nselect_bo, check nselect_bo')
        # ---------- score
        self.score = self.config.get('BO', 'score')
        if self.score in ['TS', 'EI', 'PI']:
            pass
        else:
            raise ValueError('score must be TS, EI, or PI')
        # ---------- num_rand_basis
        try:
            self.num_rand_basis = self.config.getint('BO', 'num_rand_basis')
        except configparser.NoOptionError:
            self.num_rand_basis = 0
        # ---------- cdev
        try:
            self.cdev = self.config.getfloat('BO', 'cdev')
        except configparser.NoOptionError:
            self.cdev = 0.001
        # ---------- dscrpt
        self.dscrpt = self.config.get('BO', 'dscrpt')
        # ---------- FP
        if self.dscrpt == 'FP':
            # ------ ValleOganov
            try:
                from dscribe.descriptors import ValleOganov
            except ModuleNotFoundError:
                raise ModuleNotFoundError('DScribe is required for FP. --> pip3 install dscribe')
            # ------ fp_rmax
            try:
                self.fp_rmax = self.config.getfloat('BO', 'fp_rmax')
            except configparser.NoOptionError:
                self.fp_rmax = 8.0
            if self.fp_rmax < 0:
                raise ValueError('fp_rmax < 0, check fp_rmin and fp_rmax')
            # ------ fp_npoints
            try:
                self.fp_npoints = self.config.getint('BO', 'fp_npoints')
            except configparser.NoOptionError:
                self.fp_npoints = 20
            if self.fp_npoints < 1:
                raise ValueError('fp_npoints < 1, check fp_npoints')
            # ------ fp_sigma
            try:
                self.fp_sigma = self.config.getfloat('BO', 'fp_sigma')
            except configparser.NoOptionError:
                self.fp_sigma = 0.7
            if self.fp_sigma < 0:
                raise ValueError('fp_sigma < 0, check fp_sigma')
        else:
            raise ValueError('dscrpt must be FP for now')
        # ---------- max_select_bo
        try:
            self.max_select_bo = self.config.getint('BO', 'max_select_bo')
        except configparser.NoOptionError:
            self.max_select_bo = 0
        if self.max_select_bo < 0:
            raise ValueError('max_select_bo < 0, check max_select_bo')
        # ---------- manual_select_bo
        try:
            self.manual_select_bo = self.config.get('BO', 'manual_select_bo')
            self.manual_select_bo = [int(x) for x in self.manual_select_bo.split()]
            self.manual_select_bo = tuple(self.manual_select_bo)    # list --> tuple
        except configparser.NoOptionError:
            self.manual_select_bo = None
        if self.manual_select_bo is not None:
            for cid in self.manual_select_bo:
                if not 0 <= cid < self.tot_struc:
                    raise ValueError('manual_select_bo must be non-negative int and less than tot_struc')
        # ---------- emin_bo, emax_bo
        try:
            self.emin_bo = self.config.getfloat('BO', 'emin_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.emin_bo = None
        try:
            self.emax_bo = self.config.getfloat('BO', 'emax_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.emax_bo = None
        if self.emax_bo is not None and self.emin_bo is not None:
            if self.emin_bo > self.emax_bo:
                raise ValueError('emax_bo < emin_bo, check emax_bo and emin_bo')

    def _read_laqa(self):
        # ---------- nselect_laqa
        self.nselect_laqa = self.config.getint('LAQA', 'nselect_laqa')
        if self.nselect_laqa < 1:
            raise ValueError('nselect_laqa < 1, check nselect_laqa')
        elif self.tot_struc < self.nselect_laqa:
            raise ValueError('tot_struc < nselect_laqa, check nselect_laqa')
        # ---------- wf
        try:
            self.wf = self.config.getfloat('LAQA', 'wf')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.wf = 0.1
        # ---------- ws
        try:
            self.ws = self.config.getfloat('LAQA', 'ws')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.ws = 10.0

    def _read_ea(self):
        # ---------- n_pop
        self.n_pop = self.config.getint('EA', 'n_pop')
        if self.n_pop < 1:
            raise ValueError('n_pop must be positive int')
        # ---------- n_crsov
        self.n_crsov = self.config.getint('EA', 'n_crsov')
        if self.n_crsov < 0:
            raise ValueError('n_crsov must be non-negative int')
        # ---------- n_perm
        self.n_perm = self.config.getint('EA', 'n_perm')
        if self.n_perm < 0:
            raise ValueError('n_perm must be non-negative int')
        if self.n_perm != 0 and len(self.atype) == 1:
            raise ValueError('When the number of atom type is 1, n_perm must be 0')
        # ---------- n_strain
        self.n_strain = self.config.getint('EA', 'n_strain')
        if self.n_strain < 0:
            raise ValueError('n_strain must be non-negative int')
        # ---------- n_rand
        self.n_rand = self.config.getint('EA', 'n_rand')
        if self.n_rand < 0:
            raise ValueError('n_rand must be non-negative int')
        # ---------- check n_pop for EA
        if self.algo == 'EA' and self.struc_mode not in ['mol', 'mol_bs']:
            if self.n_crsov + self.n_perm + self.n_strain + self.n_rand != self.n_pop:
                raise ValueError('n_crsov + n_perm + n_strain + n_rand must be n_pop')
        # ---------- n_elite
        self.n_elite = self.config.getint('EA', 'n_elite')
        if self.n_elite < 0:
            raise ValueError('n_elite must be non-negative int')
        # ---------- fit_reverse
        try:
            self.fit_reverse = self.config.getboolean('EA', 'fit_reverse')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.fit_reverse = False
        # ---------- n_fittest
        try:
            self.n_fittest = self.config.getint('EA', 'n_fittest')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.n_fittest = 0
        if self.n_fittest < 0:
            raise ValueError('n_fittest must be non-negative int')
        # ---------- slct_func
        self.slct_func = self.config.get('EA', 'slct_func')
        if self.slct_func not in ['TNM', 'RLT']:
            raise ValueError('slct_func must be TNM or RLT')
        # ------ TNM
        if self.slct_func == 'TNM':
            # -- t_size
            try:
                self.t_size = self.config.getint('EA', 't_size')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.t_size = 3
            if self.t_size < 2:
                raise ValueError('t_size must be greater than or equal to 2')
        # ------ RLT
        elif self.slct_func == 'RLT':
            # -- a_rlt
            try:
                self.a_rlt = self.config.getfloat('EA', 'a_rlt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.a_rlt = 10.0
            # -- b_rlt
            try:
                self.b_rlt = self.config.getfloat('EA', 'b_rlt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.b_rlt = 1.0
            # -- check
            if not 0 < self.b_rlt < self.a_rlt:
                raise ValueError('must be 0 < b_rlt < a_rlt')
        # ---------- crossover
        if self.n_crsov > 0:
            # ------ crs_lat
            try:
                self.crs_lat = self.config.get('EA', 'crs_lat')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.crs_lat = 'random'
            if self.crs_lat not in ['equal', 'random']:
                raise ValueError('crs_lat must be equal or random')
            # ------ nat_diff_tole
            try:
                self.nat_diff_tole = self.config.getint('EA', 'nat_diff_tole')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.nat_diff_tole = 4
            if self.nat_diff_tole < 0:
                raise ValueError('nat_diff_tole must be non-negative int')
        # ---------- permutation
        if self.n_perm > 0:
            # ------ ntimes
            try:
                self.ntimes = self.config.getint('EA', 'ntimes')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.ntimes = 1
            if self.ntimes < 1:
                raise ValueError('ntimes must be positive int')
        # ---------- strain
        if self.n_strain > 0:
            # ------ sigma_st
            try:
                self.sigma_st = self.config.getfloat('EA', 'sigma_st')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.sigma_st = 0.5
            if self.sigma_st <= 0:
                raise ValueError('sigma_st must be positive float')
        # ---------- other params
        # ------ maxcnt_ea
        try:
            self.maxcnt_ea = self.config.getint('EA', 'maxcnt_ea')
        except configparser.NoOptionError:
            self.maxcnt_ea = 50
        # ------ maxgen_ea
        try:
            self.maxgen_ea = self.config.getint('EA', 'maxgen_ea')
        except configparser.NoOptionError:
            self.maxgen_ea = 0
        if self.maxgen_ea < 0:
            raise ValueError('maxgen_ea must be non-negative int')
        # ------ emin_ea, emax_ea
        try:
            self.emin_ea = self.config.getfloat('EA', 'emin_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.emin_ea = None
        try:
            self.emax_ea = self.config.getfloat('EA', 'emax_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.emax_ea = None
        if self.emax_ea is not None and self.emin_ea is not None:
            if self.emin_ea > self.emax_ea:
                raise ValueError('emax_ea < emin_ea, check emax_ea and emin_ea')
        # ---------- EA-vc
        if self.algo == 'EA-vc':
            # ------ n_add
            self.n_add = self.config.getint('EA', 'n_add')
            if self.n_add < 0:
                raise ValueError('n_add must be non-negative int')
            # ------ n_elim
            self.n_elim = self.config.getint('EA', 'n_elim')
            if self.n_elim < 0:
                raise ValueError('n_elim must be non-negative int')
            # ------ n_subs
            self.n_subs = self.config.getint('EA', 'n_subs')
            if self.n_subs < 0:
                raise ValueError('n_subs must be non-negative int')
            if self.n_subs != 0:
                diff_nat = [u - l for u, l in zip(self.ul_nat, self.ll_nat)]
                non_zero_count = len([x for x in diff_nat if x != 0])
                if non_zero_count < 2:
                    raise ValueError('n_subs: at least two atomic types are variable.')
            # ------ target
            self.target = self.config.get('EA', 'target')
            #if self.target not in ['random','depop','overpop']:
            if self.target not in ['random']:
                raise ValueError('target must be random for now')
            # ------ check n_pop for EA-vc
            if self.n_crsov + self.n_perm + self.n_strain + self.n_rand \
                + self.n_add + self.n_elim + self.n_subs != self.n_pop:
                raise ValueError('n_crsov + n_perm + n_strain + n_rand'
                                    ' + n_add + n_elim + n_subs must be n_pop')
            # ------ end_point
            self.end_point = self.config.get('EA', 'end_point')
            self.end_point = tuple([float(x) for x in self.end_point.split()])
            if not len(self.end_point) == len(self.atype):
                raise ValueError('len(end_point) == len(atype), check end_point')
            # ------ show_max
            try:
                self.show_max = self.config.getfloat('EA', 'show_max')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.show_max = 0.2
            # ------ label_stable
            try:
                self.label_stable = self.config.getboolean('EA', 'label_stable')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.label_stable = True
            # ------ vmax
            try:
                self.vmax = self.config.getfloat('EA', 'vmax')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.vmax = 0.2
            # ------ bottom_margin
            try:
                self.bottom_margin = self.config.getfloat('EA', 'bottom_margin')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.bottom_margin = 0.02
        # ---------- mol or mol_bs
        if self.struc_mode in ['mol', 'mol_bs']:
            # ------ n_rotation
            self.n_rotation = self.config.getint('EA', 'n_rotation')
            if self.n_rotation < 0:
                raise ValueError('n_rotation must be non-negative int')
            # ------ check n_pop for mol or mol_bs
            if self.n_crsov + self.n_perm + self.n_strain + self.n_rand \
                + self.n_rotation != self.n_pop:
                raise ValueError('n_crsov + n_perm + n_strain + n_rand + n_rotation must be n_pop')
            # ------ mindist_mol_ea
            self.mindist_mol_ea = []
            for i in range(len(self.mol_file)):
                tmp = self.config.get('EA', f'mindist_mol_ea_{i+1}')
                tmp = tuple([float(x) for x in tmp.split()])    # str --> float
                if not len(tmp) == len(self.mol_file):
                    raise ValueError(f'not len(mindist_mol_ea_{i+1}) == len(mol_file)')
                self.mindist_mol_ea.append(tmp)
            self.mindist_mol_ea = tuple(self.mindist_mol_ea)    # list --> tuple
            # check symmetric matrix
            for i in range(len(self.mindist_mol_ea)):
                for j in range(len(self.mindist_mol_ea)):
                    if i < j:
                        if not self.mindist_mol_ea[i][j] == self.mindist_mol_ea[j][i]:
                            raise ValueError(f'mindist_mol_ea is not symmetric. ({i}, {j}) -->'
                                                f' {self.mindist_mol_ea[i][j]}, ({j}, {i}) --> {self.mindist_mol_ea[j][i]}')
            # ------ rot_max_angle
            try:
                self.rot_max_angle = self.config.getint('EA', 'rot_max_angle')
            except configparser.NoOptionError:
                self.rot_max_angle = 360
            if self.rot_max_angle <= 0:
                raise ValueError('rot_max_angle must be positive int')
            # ------ protect_mol_struc
            try:
                self.protect_mol_struc = self.config.getboolean('EA', 'protect_mol_struc')
            except configparser.NoOptionError:
                self.protect_mol_struc = False

    def _read_vasp(self):
        # ---------- kpt_flag
        self.kpt_flag = True
        # ---------- kppvol
        self.kppvol = self.config.get('VASP', 'kppvol')
        self.kppvol = tuple([int(x) for x in self.kppvol.split()])    # str --> int --> list --> tuple
        if not len(self.kppvol) == self.nstage:
            raise ValueError('not len(kppvol) == nstage, check kppvol')
        # ---------- force_gamma
        try:
            self.force_gamma = self.config.getboolean('VASP', 'force_gamma')
        except configparser.NoOptionError:
            self.force_gamma = False

    def _read_qe(self):
        # ---------- kpt_flag
        self.kpt_flag = True
        # ---------- kppvol
        self.kppvol = self.config.get('QE', 'kppvol')
        self.kppvol = tuple([int(x) for x in self.kppvol.split()])    # str --> int --> list --> tuple
        if not len(self.kppvol) == self.nstage:
            raise ValueError('not len(kppvol) == nstage, check kppvol')
        # ---------- force_gamma
        try:
            self.force_gamma = self.config.getboolean('QE', 'force_gamma')
        except configparser.NoOptionError:
            self.force_gamma = False
        # ---------- qe_infile, qe_outfile
        self.qe_infile = self.config.get('QE', 'qe_infile')
        self.qe_outfile = self.config.get('QE', 'qe_outfile')
        # ---------- pv_term
        try:
            self.pv_term = self.config.getboolean('QE', 'pv_term')
        except configparser.NoOptionError:
            self.pv_term = False
        if self.energy_step_flag and self.pv_term:
            raise ValueError('cannot parse energy_step with pv_term yet. use pv_term = False')

    def _read_omx(self):
        # ---------- kpt_flag
        self.kpt_flag = False
        # ---------- kppvol
        self.kppvol = self.config.get('OMX', 'kppvol')
        self.kppvol = tuple([int(x) for x in self.kppvol.split()])    # str --> int --> list --> tuple
        if not len(self.kppvol) == self.nstage:
            raise ValueError('not len(kppvol) == nstage, check kppvol')
        # ---------- force_gamma
        try:
            self.force_gamma = self.config.getboolean('OMX', 'force_gamma')
        except configparser.NoOptionError:
            self.force_gamma = False
        # ---------- OMX_infile, OMX_outfile
        self.OMX_infile = self.config.get('OMX', 'OMX_infile')
        self.OMX_outfile = self.config.get('OMX', 'OMX_outfile')
        # ---------- upSpin, downSpin
        self.upSpin = {}
        self.downSpin = {}
        ValenceElec = self.config.get('OMX', 'ValenceElectrons')
        ValElecIn = ValenceElec.split()
        for i in range(0, len(ValElecIn), 3):
            self.upSpin[ValElecIn[i]]   = ValElecIn[i+1]
            self.downSpin[ValElecIn[i]] = ValElecIn[i+2]

    def _read_soiap(self):
        # ---------- kpt_flag
        self.kpt_flag = False
        # ---------- force_gamma
        self.force_gamma = False
        # ---------- soiap_infile, soiap_outfile, soiap_cif
        self.soiap_infile = self.config.get('soiap', 'soiap_infile')
        self.soiap_outfile = self.config.get('soiap', 'soiap_outfile')
        self.soiap_cif = self.config.get('soiap', 'soiap_cif')

    def _read_lammps(self):
        # ---------- kpt_flag
        self.kpt_flag = False
        # ---------- force_gamma
        self.force_gamma = False
        # ---------- lammps_infile, lammps_outfile
        self.lammps_infile = self.config.get('LAMMPS', 'lammps_infile')
        self.lammps_outfile = self.config.get('LAMMPS', 'lammps_outfile')
        # ---------- lammps_potential
        try:
            self.lammps_potential = self.config.get('LAMMPS', 'lammps_potential')
            self.lammps_potential = self.lammps_potential.split()
        except configparser.NoOptionError:
            self.lammps_potential = None
        # ---------- lammps_data
        self.lammps_data = self.config.get('LAMMPS', 'lammps_data')

    def _read_ase(self):
        # ---------- kpt_flag
        self.kpt_flag = False
        # ---------- force_gamma
        self.force_gamma = False
        # ---------- ase_python
        self.ase_python = self.config.get('ASE', 'ase_python')
