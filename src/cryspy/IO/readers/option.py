import configparser
from .base import BaseReader


class OptionReader(BaseReader):
    """
    Reader for [option] section
    """
    def read(self, ht=False):
        # ---------- rslt_out
        try:
            self.rin.rslt_out = self.config.get('option', 'rslt_out')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.rslt_out = 'always'
        if self.rin.rslt_out not in ['always', 'cycle', 'off']:
            raise ValueError('rslt_out must be always, cycle, or off')

        # ---------- check_mindist_opt
        try:
            self.rin.check_mindist_opt = self.config.getboolean('option', 'check_mindist_opt')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.check_mindist_opt = True

        # ---------- stop_chkpt
        try:
            self.rin.stop_chkpt = self.config.getint('option', 'stop_chkpt')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.stop_chkpt = 0

        # ---------- load_struc_flag
        try:
            self.rin.load_struc_flag = self.config.getboolean('option', 'load_struc_flag')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.load_struc_flag = False

        # ---------- stop_next_struc
        try:
            self.rin.stop_next_struc = self.config.getboolean('option', 'stop_next_struc')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.stop_next_struc = False

        # ---------- backup_interval
        try:
            self.rin.backup_interval = self.config.getint(
                'option',
                'backup_interval',
            )
        except (configparser.NoOptionError, configparser.NoSectionError):
            if ht or self.rin.algo == 'RS':
                self.rin.backup_interval = 0
            else:
                self.rin.backup_interval = 1
        if self.rin.backup_interval < 0:
            raise ValueError(
                'backup_interval must be non-negative int'
            )

        # ---------- hull_output_interval
        try:
            self.rin.hull_output_interval = self.config.getint(
                'option',
                'hull_output_interval',
            )
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.hull_output_interval = 1
        if self.rin.hull_output_interval < 0:
            raise ValueError(
                'hull_output_interval must be non-negative int'
            )

        # ---------- recalc
        try:
            self.rin.recalc = self.config.get('option', 'recalc')
            self.rin.recalc = tuple([int(x) for x in self.rin.recalc.split()])
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.recalc = None
        if self.rin.recalc is not None:
            if self.rin.algo not in ['EA', 'EA-vc']:
                for cid in self.rin.recalc:
                    if not 0 <= cid < self.rin.tot_struc:
                        raise ValueError('recalc must be non-negative int and less than tot_struc')

        # ---------- append_struc_ea
        try:
            self.rin.append_struc_ea = self.config.getboolean('option', 'append_struc_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.append_struc_ea = False

        # ---------- energy_step_flag
        try:
            self.rin.energy_step_flag = self.config.getboolean('option', 'energy_step_flag')
            if self.rin.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('energy_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.energy_step_flag = False

        # ---------- struc_step_flag
        try:
            self.rin.struc_step_flag = self.config.getboolean('option', 'struc_step_flag')
            if self.rin.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('struc_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.struc_step_flag = False

        # ---------- force_step_flag
        try:
            self.rin.force_step_flag = self.config.getboolean('option', 'force_step_flag')
            if self.rin.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('force_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.force_step_flag = False
        if self.rin.algo == 'LAQA':
            self.rin.force_step_flag = True

        # ---------- stress_step_flag
        try:
            self.rin.stress_step_flag = self.config.getboolean('option', 'stress_step_flag')
            if self.rin.calc_code in ['LAMMPS', 'OMX', 'ASE']:
                raise ValueError('stress_step_flag: only VASP, QE, and soiap for now')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.stress_step_flag = False
        if self.rin.algo == 'LAQA':
            self.rin.stress_step_flag = True

        # ---------- seed
        try:
            self.rin.seed = self.config.getint('option', 'seed')
            if not 0 <= self.rin.seed <= 2**32 - 1:
                raise ValueError('seed must be an int between 0 and 2**32 - 1')
        except (configparser.NoOptionError, configparser.NoSectionError):
            self.rin.seed = None
