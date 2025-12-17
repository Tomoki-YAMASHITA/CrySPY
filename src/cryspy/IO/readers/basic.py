from .base import BaseReader


class BasicReader(BaseReader):
    """Reader for [basic] section"""

    def read(self):
        # ---------- algo
        self.rin.algo = self.config.get('basic', 'algo')
        if self.rin.algo not in ['RS', 'BO', 'LAQA', 'EA', 'EA-vc']:
            raise ValueError('algo must be RS, BO, LAQA, EA or EA-vc')

        # ---------- calc_code
        self.rin.calc_code = self.config.get('basic', 'calc_code')
        if self.rin.calc_code not in ['VASP', 'QE', 'soiap', 'LAMMPS', 'OMX', 'ASE']:
            raise ValueError('calc_code must be VASP, QE, OMX, soiap, LAMMPS, or ASE')
        if self.rin.algo == 'LAQA':
            if self.rin.calc_code not in ['VASP', 'QE', 'soiap']:
                raise ValueError('LAQA: only VASP, QE, or soiap for now')

        # ---------- tot_struc
        if self.rin.algo not in ['EA', 'EA-vc']:
            self.rin.tot_struc = self.config.getint('basic', 'tot_struc')
            if self.rin.tot_struc < 1:
                raise ValueError('tot_struc < 1, check tot_struc')

        # ---------- nstage
        self.rin.nstage = self.config.getint('basic', 'nstage')
        if self.rin.nstage < 1:
            raise ValueError('nstage < 1, check nstage')
        if self.rin.algo == 'LAQA':
            if not self.rin.nstage == 1:
                raise ValueError('nstage must be 1 in LAQA')

        # ---------- njob
        self.rin.njob = self.config.getint('basic', 'njob')
        if self.rin.njob < 1:
            raise ValueError('njob < 1, check njob')

        # ---------- jobcmd
        self.rin.jobcmd = self.config.get('basic', 'jobcmd')

        # ---------- jobfile
        self.rin.jobfile = self.config.get('basic', 'jobfile')