import configparser
from .base import BaseReader


class SoiapReader(BaseReader):
    """
    Reader for [soiap] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = False

        # ---------- force_gamma
        self.rin.force_gamma = False

        # ---------- soiap_infile, soiap_outfile, soiap_cif
        self.rin.soiap_infile = self.config.get('soiap', 'soiap_infile')
        self.rin.soiap_outfile = self.config.get('soiap', 'soiap_outfile')
        self.rin.soiap_cif = self.config.get('soiap', 'soiap_cif')
