import configparser
from .base import BaseReader


class LAQAReader(BaseReader):
    """
    Reader for [LAQA] section
    """
    def read(self):
        # ---------- nselect_laqa
        self.rin.nselect_laqa = self.config.getint('LAQA', 'nselect_laqa')
        if self.rin.nselect_laqa < 1:
            raise ValueError('nselect_laqa < 1, check nselect_laqa')
        elif self.rin.tot_struc < self.rin.nselect_laqa:
            raise ValueError('tot_struc < nselect_laqa, check nselect_laqa')

        # ---------- wf
        try:
            self.rin.wf = self.config.getfloat('LAQA', 'wf')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.wf = 0.1

        # ---------- ws
        try:
            self.rin.ws = self.config.getfloat('LAQA', 'ws')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.ws = 10.0