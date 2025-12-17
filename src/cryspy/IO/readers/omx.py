import configparser
from .base import BaseReader


class OMXReader(BaseReader):
    """
    Reader for [OMX] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = False

        # ---------- kppvol
        kppvol_str = self.config.get('OMX', 'kppvol')
        self.rin.kppvol = tuple(int(x) for x in kppvol_str.split())
        if len(self.rin.kppvol) != self.rin.nstage:
            raise ValueError('not len(kppvol) == nstage, check OMX.kppvol')

        # ---------- force_gamma
        try:
            self.rin.force_gamma = self.config.getboolean('OMX', 'force_gamma')
        except configparser.NoOptionError:
            self.rin.force_gamma = False

        # ---------- OMX_infile, OMX_outfile
        self.rin.OMX_infile = self.config.get('OMX', 'OMX_infile')
        self.rin.OMX_outfile = self.config.get('OMX', 'OMX_outfile')

        # ---------- upSpin, downSpin
        self.rin.upSpin = {}
        self.rin.downSpin = {}
        valence = self.config.get('OMX', 'ValenceElectrons')
        tokens = valence.split()
        if len(tokens) % 3 != 0:
            raise ValueError('ValenceElectrons must be triplets: "Elem up down"')
        for i in range(0, len(tokens), 3):
            elem = tokens[i]
            self.rin.upSpin[elem] = tokens[i + 1]
            self.rin.downSpin[elem] = tokens[i + 2]