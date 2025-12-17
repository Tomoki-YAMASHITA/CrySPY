import configparser
from .base import BaseReader


class QEReader(BaseReader):
    """
    Reader for [QE] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = True

        # ---------- kppvol (required)
        try:
            kppvol_str = self.config.get('QE', 'kppvol')
        except (configparser.NoOptionError, configparser.NoSectionError) as e:
            raise ValueError('QE.kppvol is required') from e
        self.rin.kppvol = tuple(int(x) for x in kppvol_str.split())
        if len(self.rin.kppvol) != self.rin.nstage:
            raise ValueError('not len(kppvol) == nstage, check QE.kppvol')

        # ---------- force_gamma
        try:
            self.rin.force_gamma = self.config.getboolean('QE', 'force_gamma')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.force_gamma = False

        # ---------- qe_infile, qe_outfile (required)
        try:
            self.rin.qe_infile = self.config.get('QE', 'qe_infile')
            self.rin.qe_outfile = self.config.get('QE', 'qe_outfile')
        except (configparser.NoOptionError, configparser.NoSectionError) as e:
            raise ValueError('QE.qe_infile and QE.qe_outfile are required') from e

        # ---------- pv_term
        try:
            self.rin.pv_term = self.config.getboolean('QE', 'pv_term')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.pv_term = False
        if self.rin.energy_step_flag and self.rin.pv_term:
            raise ValueError('cannot parse energy_step with pv_term yet. use pv_term = False')