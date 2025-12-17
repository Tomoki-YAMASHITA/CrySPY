import configparser
from .base import BaseReader


class VASPReader(BaseReader):
    """
    Reader for [VASP] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = True

        # ---------- kppvol
        try:
            kppvol_str = self.config.get('VASP', 'kppvol')
        except (configparser.NoOptionError, configparser.NoSectionError):
            raise ValueError('VASP.kppvol is required')
        self.rin.kppvol = tuple(int(x) for x in kppvol_str.split())
        if len(self.rin.kppvol) != self.rin.nstage:
            raise ValueError('not len(kppvol) == nstage, check VASP.kppvol')

        # ---------- force_gamma
        try:
            self.rin.force_gamma = self.config.getboolean('VASP', 'force_gamma')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.force_gamma = False

        # ---------- EA-vc specific
        if self.rin.algo in ['EA-vc']:
            # ------ vasp_MAGMOM
            try:
                s = self.config.get('VASP', 'vasp_MAGMOM')
                self.rin.vasp_MAGMOM = tuple(float(a) for a in s.split())
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.vasp_MAGMOM = None
            # ------ vasp_LDAUL
            try:
                s = self.config.get('VASP', 'vasp_LDAUL')
                self.rin.vasp_LDAUL = tuple(int(a) for a in s.split())
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.vasp_LDAUL = None
            # ------ vasp_LDAUU
            try:
                s = self.config.get('VASP', 'vasp_LDAUU')
                self.rin.vasp_LDAUU = tuple(float(a) for a in s.split())
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.vasp_LDAUU = None
            # ------ vasp_LDAUJ
            try:
                s = self.config.get('VASP', 'vasp_LDAUJ')
                self.rin.vasp_LDAUJ = tuple(float(a) for a in s.split())
            except (configparser.NoOptionError, configparser.NoSectionError):
                self.rin.vasp_LDAUJ = None