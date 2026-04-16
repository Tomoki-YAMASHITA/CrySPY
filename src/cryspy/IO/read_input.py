import configparser
from dataclasses import dataclass, field
from logging import getLogger
import os


logger = getLogger('cryspy')


@dataclass
class ReadInput:
    '''
    Read cryspy.in
        Here, ignore type hint for None
    '''
    # ---------- basic section
    algo: str = field(default=None)
    calc_code: str = field(default=None)
    tot_struc: int = field(default=None)
    nstage: int = field(default=None)
    njob: int = field(default=None)
    jobcmd: str = field(default=None)
    jobfile: str = field(default=None)

    # ---------- structure section
    struc_mode: str = field(default=None)
    atype: tuple = field(default=None)
    nat: tuple = field(default=None)
    mindist: tuple = field(default=None)
    mindist_factor: float = field(default=None)
    vol_factor: float = field(default=None)
    vol_mu: float = field(default=None)
    vol_sigma: float = field(default=None)
    symprec: float = field(default=None)
    spgnum: str = field(default=None)    # ignore type hint here. 'all', tuple or 0
    use_find_wy: bool = field(default=None)
    # ------ EA-vc
    ll_nat: tuple = field(default=None)
    ul_nat: tuple = field(default=None)
    charge: tuple = field(default=None)
    min_comp: tuple = field(default=None)
    max_comp: tuple = field(default=None)
    # ------ mol or mol_bs
    mol_file: tuple = field(default=None)
    nmol: tuple = field(default=None)
    timeout_mol: float = field(default=None)
    rot_mol: str = field(default=None)
    nrot: int = field(default=None)
    mindist_mol_bs: tuple = field(default=None)
    mindist_mol_bs_factor: float = field(default=None)
    # ------ find_wy or spgnum = 0
    fwpath: str = field(default=None)
    minlen: float = field(default=None)
    maxlen: float = field(default=None)
    dangle: float = field(default=None)
    maxcnt: int = field(default=None)

    # ---------- visual section
    ymax: float = field(default=None)
    markersize: int = field(default=None)
    fig_format: str = field(default=None)
    # ------ for EA
    plot_min_gen: int = field(default=None)
    plot_max_gen: int = field(default=None)
    ref_gen: int = field(default=None)
    # ------ for EA-vc
    show_max: float = field(default=None)
    label_stable: bool = field(default=None)
    vmax: float = field(default=None)
    bottom_margin: float = field(default=None)
    axis_order: str = field(default=None)

    # ---------- option section
    check_mindist_opt: bool = field(default=None)
    stop_chkpt: int = field(default=None)
    load_struc_flag: bool = field(default=None)
    stop_next_struc: bool = field(default=None)
    recalc: tuple = field(default=None)
    append_struc_ea: bool = field(default=None)
    energy_step_flag: bool = field(default=None)
    struc_step_flag: bool = field(default=None)
    force_step_flag: bool = field(default=None)
    stress_step_flag: bool = field(default=None)
    seed: int = field(default=None)

    # ---------- BO section
    nselect_bo: int = field(default=None)
    score: str = field(default=None)
    num_rand_basis: int = field(default=None)
    cdev: float = field(default=None)
    dscrpt: str = field(default=None)
    fp_rmin: float = field(default=None)
    fp_rmax: float = field(default=None)
    fp_npoints: int = field(default=None)
    fp_sigma: float = field(default=None)
    max_select_bo: int = field(default=None)
    manual_select_bo: tuple = field(default=None)
    emin_bo: float = field(default=None)
    emax_bo: float = field(default=None)

    # ---------- LAQA section
    nselect_laqa: int = field(default=None)
    wf: float = field(default=None)
    ws: float = field(default=None)

    # ---------- EA section
    n_pop: int = field(default=None)
    n_crsov: int = field(default=None)
    n_perm: int = field(default=None)
    n_strain: int = field(default=None)
    n_rand: int = field(default=None)
    n_elite: int = field(default=None)
    fit_reverse: bool = field(default=None)
    n_fittest: int = field(default=None)
    slct_func: str = field(default=None)
    t_size: int = field(default=None)
    a_rlt: float = field(default=None)
    b_rlt: float = field(default=None)
    crs_lat: str = field(default=None)
    nat_diff_tole: int = field(default=None)
    ntimes: int = field(default=None)
    sigma_st: float = field(default=None)
    maxcnt_ea: int = field(default=None)
    maxgen_ea: int = field(default=None)
    emin_ea: float = field(default=None)
    emax_ea: float = field(default=None)
    n_add: int = field(default=None)
    add_max: int = field(default=None)
    n_elim: int = field(default=None)
    elim_max: int = field(default=None)
    n_subs: int = field(default=None)
    subs_max: int = field(default=None)
    target: str = field(default=None)
    end_point: tuple = field(default=None)
    n_rotation: int = field(default=None)          # not implemented yet, for EA mol
    mindist_mol_ea: tuple = field(default=None)    # not implemented yet, for EA mol
    rot_max_angle: float = field(default=None)     # not implemented yet, for EA mol
    protect_mol_struc: bool = field(default=None)  # not implemented yet, for EA mol

    # ---------- common in VASP, QE, OMX
    kpt_flag: bool = field(default=None)
    kppvol: int = field(default=None)
    force_gamma: bool = field(default=None)

    # ---------- VASP section
    vasp_MAGMOM: tuple = field(default=None)
    vasp_LDAUL: tuple = field(default=None)
    vasp_LDAUU: tuple = field(default=None)
    vasp_LDAUJ: tuple = field(default=None)

    # ---------- QE section
    qe_infile: str = field(default=None)
    qe_outfile: str = field(default=None)
    pv_term: bool = field(default=None)

    # ---------- OMX section
    OMX_infile: str = field(default=None)
    OMX_outfile: str = field(default=None)
    upSpin: dict = field(default=None)
    downSpin: dict = field(default=None)

    # ---------- soiap section
    soiap_infile: str = field(default=None)
    soiap_outfile: str = field(default=None)
    soiap_cif: str = field(default=None)

    # ---------- lammps section
    lammps_infile: str = field(default=None)
    lammps_outfile: str = field(default=None)
    lammps_potential: str = field(default=None)
    lammps_data: str = field(default=None)

    # ---------- ASE section
    ase_python: str = field(default=None)

    def __post_init__(self):
        # ---------- import readers
        from .readers.basic import BasicReader
        from .readers.structure import StructureReader
        from .readers.visual import VisualReader
        from .readers.option import OptionReader
        from .readers.bo import BOReader
        from .readers.laqa import LAQAReader
        from .readers.ea import EAReader
        from .readers.vasp import VASPReader
        from .readers.qe import QEReader
        from .readers.omx import OMXReader
        from .readers.soiap import SoiapReader
        from .readers.lammps import LAMMPSReader
        from .readers.ase import ASEReader

        # ---------- read cryspy.in
        filename = 'cryspy.in'
        if not os.path.isfile(filename):
            raise FileNotFoundError(f'{filename} not found in {os.getcwd()}')
        config = configparser.ConfigParser()
        config.read(filename)

        # ---------- read sections
        BasicReader(config, self).read()
        StructureReader(config, self).read()
        VisualReader(config, self).read()
        OptionReader(config, self).read()    # do this prior to EA for append_struc_ea
        if self.algo == 'BO':
            BOReader(config, self).read()
        if self.algo == 'LAQA':
            LAQAReader(config, self).read()
        if self.algo in ['EA', 'EA-vc'] or self.append_struc_ea:
            EAReader(config, self).read()
        if self.calc_code == 'VASP':
            VASPReader(config, self).read()
        if self.calc_code == 'QE':
            QEReader(config, self).read()
        if self.calc_code == 'OMX':
            OMXReader(config, self).read()
        if self.calc_code == 'soiap':
            SoiapReader(config, self).read()
        if self.calc_code == 'LAMMPS':
            LAMMPSReader(config, self).read()
        if self.calc_code == 'ASE':
            ASEReader(config, self).read()
