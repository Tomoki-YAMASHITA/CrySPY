import configparser
from .base import BaseReader


class LAMMPSReader(BaseReader):
    """
    Reader for [LAMMPS] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = False

        # ---------- force_gamma
        self.rin.force_gamma = False

        # ---------- lammps_infile, lammps_outfile
        self.rin.lammps_infile = self.config.get('LAMMPS', 'lammps_infile')
        self.rin.lammps_outfile = self.config.get('LAMMPS', 'lammps_outfile')

        # ---------- lammps_potential
        try:
            self.rin.lammps_potential = self.config.get('LAMMPS', 'lammps_potential')
            self.rin.lammps_potential = self.rin.lammps_potential.split()
        except configparser.NoOptionError:
            self.rin.lammps_potential = None

        # ---------- lammps_data
        self.rin.lammps_data = self.config.get('LAMMPS', 'lammps_data')