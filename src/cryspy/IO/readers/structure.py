import configparser
import re
from .base import BaseReader

# ---------- import later
#from ...util.struc_util import get_feasible_composition


class StructureReader(BaseReader):
    """
    Reader for [structure] section
    """
    def read(self):
        # ---------- struc_mode
        try:
            self.rin.struc_mode = self.config.get('structure', 'struc_mode')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.struc_mode = 'crystal'
        if self.rin.struc_mode not in ['crystal', 'mol', 'mol_bs']:
            raise ValueError('struc_mode must be crystal, mol, or  mol_bs')
        if self.rin.algo in ['EA', 'EA-vc'] and self.rin.struc_mode in ['mol', 'mol_bs']:
            raise ValueError('currently, molecular structure is not allowed in EA or EA-vc')

        # ---------- atype
        self.rin.atype = self.config.get('structure', 'atype')
        self.rin.atype = tuple([a for a in self.rin.atype.split()])    # str --> list --> tuple
        if self.rin.algo == 'EA-vc':
            if len(self.rin.atype) < 2:
                raise ValueError('EA-vc: atype must have at least 2 elements')

        # ---------- nat
        if not self.rin.algo == 'EA-vc':
            self.rin.nat = self.config.get('structure', 'nat')
            self.rin.nat = tuple([int(x) for x in self.rin.nat.split()])    # str --> int --> list --> tuple
            if not len(self.rin.nat) == len(self.rin.atype):
                raise ValueError('not len(nat) == len(atype), check atype and nat')

        # ---------- mindist
        try:
            self.rin.mindist = []
            for i in range(len(self.rin.atype)):
                tmp = self.config.get('structure', f'mindist_{i+1}')
                tmp = tuple([float(x) for x in tmp.split()])    # str --> float --> list --> tuple
                if not len(tmp) == len(self.rin.atype):
                    raise ValueError(f'not len(mindist_{i+1}) == len(atype)')
                self.rin.mindist.append(tmp)
            self.rin.mindist = tuple(self.rin.mindist)    # list --> tuple
            # ------ check symmetric matrix
            for i in range(len(self.rin.mindist)):
                for j in range(len(self.rin.mindist)):
                    if i < j:
                        if not self.rin.mindist[i][j] == self.rin.mindist[j][i]:
                            raise ValueError(f'mindist is not symmetric. ({i}, {j}) -->'
                                            f' {self.rin.mindist[i][j]}, ({j}, {i}) --> {self.rin.mindist[j][i]}')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.mindist = None

        # ---------- mindist_factor
        try:
            self.rin.mindist_factor = self.config.getfloat('structure', 'mindist_factor')
            if self.rin.mindist_factor <= 0.0:
                raise ValueError('mindist_factor must be positive')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.mindist_factor = 1.0

        # ---------- vol_factor
        try:
            self.rin.vol_factor = self.config.getfloat('structure', 'vol_factor')
            if self.rin.vol_factor <= 0.0:
                raise ValueError('vol_factor must be positive')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.vol_factor = 1.1

        # ---------- vol_mu
        try:
            self.rin.vol_mu = self.config.getfloat('structure', 'vol_mu')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.vol_mu = None
        if self.rin.vol_mu is not None:
            if self.rin.vol_mu <= 0.0:
                raise ValueError('vol_mu must be positive float')

        # ---------- vol_sigma
        try:
            self.rin.vol_sigma = self.config.getfloat('structure', 'vol_sigma')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.vol_sigma = None
        if self.rin.vol_mu is not None and self.rin.vol_sigma is None:
            raise ValueError('vol_mu is set, but vol_sigma is not set')
        if self.rin.vol_sigma is not None:
            if self.rin.vol_sigma < 0.0:
                raise ValueError('vol_sigma must not be negative')

        # ---------- symprec
        try:
            self.rin.symprec = self.config.getfloat('structure', 'symprec')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.symprec = 0.01
        if self.rin.symprec < 0.0:
            raise ValueError('symprec must be positive float')

        # ---------- spgnum
        try:
            self.rin.spgnum = self.config.get('structure', 'spgnum')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.spgnum = 'all'
        if self.rin.spgnum == '0':    # str here
            if self.rin.struc_mode in ['mol', 'mol_bs']:
                raise ValueError('spgnum = 0 is not allow when struc_mode is mol or mol_bs')
            self.rin.spgnum = 0    # int here
        elif self.rin.spgnum == 'all':
            pass
        else:
            self._spgtuple()    # self.rin.spgnum: e.g. 3 15-18 100 --> (3, 15, 16, 17, 18, 100)

        # ---------- use_find_wy
        try:
            self.rin.use_find_wy = self.config.getboolean('structure', 'use_find_wy')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.use_find_wy = False
        if self.rin.use_find_wy:
            if not self.rin.struc_mode == 'crystal':
                raise ValueError('find_wy can be use if struc_mode is crystal')

        # ---------- EA-vc
        if self.rin.algo == 'EA-vc':
            # ------ ll_nat, ul_nat
            self.rin.ll_nat = self.config.get('structure', 'll_nat')
            self.rin.ll_nat = tuple([int(x) for x in self.rin.ll_nat.split()])    # str --> int --> list --> tuple
            self.rin.ul_nat = self.config.get('structure', 'ul_nat')
            self.rin.ul_nat = tuple([int(x) for x in self.rin.ul_nat.split()])    # str --> int --> list --> tuple
            if not len(self.rin.atype) == len(self.rin.ll_nat) == len(self.rin.ul_nat):
                raise ValueError('not len(atype) == len(ll_nat) == len(ul_nat), check ll_nat and ul_nat')
            for i in range(len(self.rin.ll_nat)):
                if not 0 <= self.rin.ll_nat[i] <= self.rin.ul_nat[i]:
                    raise ValueError(f'not 0 <= ll_nat[{i}] <= ul_nat[{i}], check ll_nat and ul_nat')
            # ------ charge
            try:
                self.rin.charge = self.config.get('structure', 'charge')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.charge = None
            if self.rin.charge is not None:
                self.rin.charge = self._parse_charge(self.rin.charge)    # str --> tuple
                if len(self.rin.charge) != len(self.rin.atype):
                    raise ValueError('not len(charge) == len(atype), check charge')
                # check that there is at least one positive and one negative value
                has_pos = False
                has_neg = False
                for c in self.rin.charge:
                    if isinstance(c, int):    # monovalent
                        if c > 0:
                            has_pos = True
                        elif c < 0:
                            has_neg = True
                    else:    # multivalent
                        for v in c:
                            if v > 0:
                                has_pos = True
                            elif v < 0:
                                has_neg = True
                    if has_pos and has_neg:
                        break
                if not (has_pos and has_neg):
                    raise ValueError(
                        'charge must contain at least one positive and one negative value '
                        '(zero is allowed).'
                    )
            # ------ min_comp, max_comp
            try:
                self.rin.min_comp = self.config.get('structure', 'min_comp')
                self.rin.min_comp = tuple([float(x) for x in self.rin.min_comp.split()])
                if len(self.rin.min_comp) != len(self.rin.atype):
                    raise ValueError('not len(min_comp) == len(atype), check min_comp')
                if any((c < 0.0) or (c > 1.0) for c in self.rin.min_comp):
                    raise ValueError('min_comp must be between 0.0 and 1.0')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.min_comp = None
            try:
                self.rin.max_comp = self.config.get('structure', 'max_comp')
                self.rin.max_comp = tuple([float(x) for x in self.rin.max_comp.split()])
                if len(self.rin.max_comp) != len(self.rin.atype):
                    raise ValueError('not len(max_comp) == len(atype), check max_comp')
                if any((c < 0.0) or (c > 1.0) for c in self.rin.max_comp):
                    raise ValueError('max_comp must be between 0.0 and 1.0')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.max_comp = None
            # -- check min_comp and max_comp
            k = len(self.rin.atype)
            if (self.rin.min_comp is None) and (self.rin.max_comp is None):
                pass
            else:
                # allow one-sided input by filling the other side
                if self.rin.min_comp is None:
                    self.rin.min_comp = tuple(0.0 for _ in range(k))
                if self.rin.max_comp is None:
                    self.rin.max_comp = tuple(1.0 for _ in range(k))
                # check feasibility
                from ...util.struc_util import get_feasible_composition, precompute_feasible_N
                feasible_comp = get_feasible_composition(self.rin.min_comp, self.rin.max_comp)
                if feasible_comp is None:
                    raise ValueError('No feasible composition exists for the given min_comp and max_comp')
                feasible_N = precompute_feasible_N(
                    self.rin.ll_nat,
                    self.rin.ul_nat,
                    feasible_comp,
                )
                if len(feasible_N) == 0:
                    raise ValueError(
                        'No feasible nat exists for the given min_comp, max_comp, ll_nat, and ul_nat'
                    )

        # ---------- mol or mol_bs
        if self.rin.struc_mode in ['mol', 'mol_bs']:
            # ------ mol_file
            self.rin.mol_file = self.config.get('structure', 'mol_file')
            self.rin.mol_file = tuple([a for a in self.rin.mol_file.split()])    # str --> list --> tuple
            # ------ nmol
            self.rin.nmol = self.config.get('structure', 'nmol')
            self.rin.nmol = tuple([int(x) for x in self.rin.nmol.split()])    # str --> int --> list --> tuple
            if not len(self.rin.mol_file) == len(self.rin.nmol):
                raise ValueError('not len(mol_file) == len(nmol)')
            # ------ timeout_mol
            try:
                self.rin.timeout_mol = self.config.getfloat('structure', 'timeout_mol')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.timeout_mol = None
            if self.rin.timeout_mol is not None:
                if self.rin.timeout_mol <= 0:
                    raise ValueError('timeout_mol must be None or positive')

        # ---------- mol_bs
        if self.rin.struc_mode == 'mol_bs':
            # ------ rot_mol
            try:
                self.rin.rot_mol = self.config.get('structure', 'rot_mol')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.rot_mol = 'random_wyckoff'
            if self.rin.rot_mol not in ['random', 'random_mol', 'random_wyckoff']:
                raise ValueError('rot_mol must be random, random_mol, or random_wyckoff')
            # ------ nrot
            try:
                self.rin.nrot = self.config.getint('structure', 'nrot')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.nrot = 20
            if self.rin.nrot < 1:
                raise ValueError('nrot < 1, check nrot')
            # ------ mindist_mol_bs
            try:
                self.rin.mindist_mol_bs = []
                for i in range(len(self.rin.mol_file)):
                    tmp = self.config.get('structure', f'mindist_mol_bs_{i+1}')
                    tmp = tuple([float(x) for x in tmp.split()])    # str --> float --> list --> tuple
                    if not len(tmp) == len(self.rin.mol_file):
                        raise ValueError(f'not len(mindist_mol_bs_{i+1}) == len(mol_file)')
                    self.rin.mindist_mol_bs.append(tmp)
                self.rin.mindist_mol_bs = tuple(self.rin.mindist_mol_bs)    # list --> tuple
                # -- check symmetric matrix
                for i in range(len(self.rin.mindist_mol_bs)):
                    for j in range(len(self.rin.mindist_mol_bs)):
                        if i < j:
                            if not self.rin.mindist_mol_bs[i][j] == self.rin.mindist_mol_bs[j][i]:
                                raise ValueError(f'mindist_mol_bs is not symmetric. ({i}, {j}) -->'
                                                f' {self.rin.mindist_mol_bs[i][j]}, ({j}, {i}) --> {self.rin.mindist_mol_bs[j][i]}')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.mindist_mol_bs = None
            # ------ mindist_mol_bs_factor
            try:
                self.rin.mindist_mol_bs_factor = self.config.getfloat('structure', 'mindist_mol_bs_factor')
                if self.rin.mindist_mol_bs_factor <= 0.0:
                    raise ValueError('mindist_mol_bs_factor must be positive')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.mindist_mol_bs_factor = 1.0

        # ---------- find_wy or spgnum = 0
        if self.rin.spgnum == 0 or self.rin.use_find_wy:
            # ------ fwpath
            try:
                self.rin.fwpath = self.config.get('structure', 'fwpath')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.fwpath = None
            # ------ minlen
            self.rin.minlen = self.config.getfloat('structure', 'minlen')
            if self.rin.minlen <= 0.0:
                raise ValueError('minlen must be positive')
            # ------ maxlen
            self.rin.maxlen = self.config.getfloat('structure', 'maxlen')
            if self.rin.minlen > self.rin.maxlen:
                raise ValueError('minlen > maxlen')
            # ------ dangle
            self.rin.dangle = self.config.getfloat('structure', 'dangle')
            if self.rin.dangle <= 0.0:
                raise ValueError('dangle must be positive')
            # ------ maxcnt
            try:
                self.rin.maxcnt = self.config.getint('structure', 'maxcnt')
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.maxcnt = 50
            if self.rin.maxcnt < 1:
                raise ValueError('maxcnt must be positive int')

    def _spgtuple(self):
        tmpspg = []
        for c in self.rin.spgnum.split():
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
        self.rin.spgnum = tuple(tmpspg)

    def _parse_charge(self, s: str):
        """
        Parse charge string:
            "1 2 -1"        -> (1, 2, -1)
            "(2 3) -2"      -> ((2,3), -2)
            "[2, 3] -2"     -> ((2,3), -2)
            "(2) 5"         -> ERROR
        """
        toks = re.findall(r'\([^\)]*\)|\[[^\]]*\]|\S+', s)
        out = []
        for tok in toks:
            tok = tok.strip()
            # tuple case
            if tok.startswith("(") or tok.startswith("["):
                inner = tok[1:-1]
                vals = re.split(r'[,\s]+', inner)
                vals = [v for v in vals if v]   # remove empty
                if len(vals) < 2:
                    raise ValueError(
                        f"Invalid multivalence format '{tok}': "
                        "use parentheses only for 2 or more values"
                    )
                out.append(tuple(int(v) for v in vals))
            else:
                # integer case
                out.append(int(tok))
        return tuple(out)
