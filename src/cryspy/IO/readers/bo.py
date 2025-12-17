import configparser
from logging import getLogger
from .base import BaseReader


logger = getLogger('cryspy')


class BOReader(BaseReader):
    """
    Reader for [BO] section
    """
    def read(self):
        # ---------- check physbo
        try:
            import physbo
        except ModuleNotFoundError:
            raise ModuleNotFoundError('PHYSBO is required for BO. --> pip3 install physbo')

        # ---------- nselect_bo
        self.rin.nselect_bo = self.config.getint('BO', 'nselect_bo')
        if self.rin.nselect_bo < 1:
            raise ValueError('nselect_bo < 1, check nselect_bo')
        elif self.rin.tot_struc < self.rin.nselect_bo:
            raise ValueError('tot_struc < nselect_bo, check nselect_bo')

        # ---------- score
        self.rin.score = self.config.get('BO', 'score')
        if self.rin.score in ['TS', 'EI', 'PI']:
            pass
        else:
            raise ValueError('score must be TS, EI, or PI')

        # ---------- num_rand_basis
        try:
            self.rin.num_rand_basis = self.config.getint('BO', 'num_rand_basis')
        except configparser.NoOptionError:
            self.rin.num_rand_basis = 0

        # ---------- cdev
        try:
            self.rin.cdev = self.config.getfloat('BO', 'cdev')
        except configparser.NoOptionError:
            self.rin.cdev = 0.001

        # ---------- dscrpt
        self.rin.dscrpt = self.config.get('BO', 'dscrpt')

        # ---------- FP
        if self.rin.dscrpt == 'FP':
            # ------ check dscribe
            try:
                from dscribe.descriptors import ValleOganov
            except ModuleNotFoundError:
                raise ModuleNotFoundError('DScribe is required for FP. --> pip3 install dscribe')
            # ------ fp_rmax
            try:
                self.rin.fp_rmax = self.config.getfloat('BO', 'fp_rmax')
            except configparser.NoOptionError:
                self.rin.fp_rmax = 8.0
            if self.rin.fp_rmax < 0:
                raise ValueError('fp_rmax < 0, check fp_rmin and fp_rmax')
            # ------ fp_npoints
            try:
                self.rin.fp_npoints = self.config.getint('BO', 'fp_npoints')
            except configparser.NoOptionError:
                self.rin.fp_npoints = 20
            if self.rin.fp_npoints < 1:
                raise ValueError('fp_npoints < 1, check fp_npoints')
            # ------ fp_sigma
            try:
                self.rin.fp_sigma = self.config.getfloat('BO', 'fp_sigma')
            except configparser.NoOptionError:
                self.rin.fp_sigma = 0.7
            if self.rin.fp_sigma < 0:
                raise ValueError('fp_sigma < 0, check fp_sigma')
        else:
            raise ValueError('dscrpt must be FP for now')

        # ---------- max_select_bo
        try:
            self.rin.max_select_bo = self.config.getint('BO', 'max_select_bo')
        except configparser.NoOptionError:
            self.rin.max_select_bo = 0
        if self.rin.max_select_bo < 0:
            raise ValueError('max_select_bo < 0, check max_select_bo')

        # ---------- manual_select_bo
        try:
            manual_select_bo_str = self.config.get('BO', 'manual_select_bo')
            self.rin.manual_select_bo = tuple([int(x) for x in manual_select_bo_str.split()])
        except configparser.NoOptionError:
            self.rin.manual_select_bo = None
        if self.rin.manual_select_bo is not None:
            for cid in self.rin.manual_select_bo:
                if not 0 <= cid < self.rin.tot_struc:
                    raise ValueError('manual_select_bo must be non-negative int and less than tot_struc')

        # ---------- emin_bo, emax_bo
        try:
            self.rin.emin_bo = self.config.getfloat('BO', 'emin_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.emin_bo = None
        try:
            self.rin.emax_bo = self.config.getfloat('BO', 'emax_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.emax_bo = None
        if self.rin.emax_bo is not None and self.rin.emin_bo is not None:
            if self.rin.emin_bo > self.rin.emax_bo:
                raise ValueError('emax_bo < emin_bo, check emax_bo and emin_bo')