from .base import BaseReader


class ASEReader(BaseReader):
    """
    Reader for [ASE] section
    """
    def read(self):
        # ---------- kpt_flag
        self.rin.kpt_flag = False

        # ---------- force_gamma
        self.rin.force_gamma = False

        # ---------- ase_python
        self.rin.ase_python = self.config.get('ASE', 'ase_python')