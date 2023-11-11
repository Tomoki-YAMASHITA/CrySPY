'''
Read input from cryspy.in
'''

import configparser
from logging import getLogger
import os

from . import io_stat
from ..util import utility


logger = getLogger('cryspy')

def readin():
    # ---------- read cryspy.in
    if not os.path.isfile('cryspy.in'):
        logger.error('Could not find cryspy.in file')
        raise SystemExit(1)
    config = configparser.ConfigParser()
    config.read('cryspy.in')

    # ---------- basic
    # ------ global declaration
    global algo, calc_code, tot_struc
    global nstage, njob, jobcmd, jobfile
    # ------ read intput variables
    calc_code = config.get('basic', 'calc_code')
    if calc_code not in ['VASP', 'QE', 'soiap', 'LAMMPS', 'OMX', 'ASE', 'ext']:
        logger.error('calc_code must be VASP, QE, OMX, soiap, LAMMPS, ASE, or ext')
        raise SystemExit(1)
    algo = config.get('basic', 'algo')
    if algo not in ['RS', 'BO', 'LAQA', 'EA', 'EA-vc']:
        logger.error('algo must be RS, BO, LAQA, EA or EA-vc')
        raise SystemExit(1)
    if algo == 'LAQA':
        if calc_code not in ['VASP', 'QE', 'soiap']:
            logger.error('LAQA: only VASP, QE, and soiap for now')
            raise SystemExit(1)
    tot_struc = config.getint('basic', 'tot_struc')
    if tot_struc <= 0:
        logger.error('tot_struc <= 0, check tot_struc')
        raise SystemExit(1)
    if not calc_code == 'ext':
        nstage = config.getint('basic', 'nstage')
        if nstage <= 0:
            logger.error('nstage <= 0, check nstage')
            raise SystemExit(1)
        if algo == 'LAQA':
            if not nstage == 1:
                logger.error('nstage shoud be 1 in LAQA')
                raise SystemExit(1)
        njob = config.getint('basic', 'njob')
        if njob <= 0:
            logger.error('njob <= 0, check njob')
            raise SystemExit(1)
        jobcmd = config.get('basic', 'jobcmd')
        jobfile = config.get('basic', 'jobfile')

    # ---------- structure
    # ------ global declaration
    global struc_mode, natot, atype, nat, ll_nat, ul_nat
    global mol_file, nmol, timeout_mol, rot_mol, nrot
    global vol_factor, vol_mu, vol_sigma
    global mindist, mindist_factor, mindist_mol_bs, mindist_mol_bs_factor
    global maxcnt, symprec, spgnum, use_find_wy, fwpath
    global minlen, maxlen, dangle

    # ------ read intput variables
    try:
        struc_mode = config.get('structure', 'struc_mode')
    except (configparser.NoOptionError, configparser.NoSectionError):
        struc_mode = 'crystal'
    if struc_mode not in ['crystal', 'mol', 'mol_bs', 'host']:
        logger.error('struc_mode is wrong')
        raise SystemExit(1)
    if algo in ['EA', 'EA-vc'] and struc_mode in ['mol', 'mol_bs']:
        if calc_code == 'ext':
            logger.error('EA, mol or mol_bs, ext, not yet')
            raise SystemExit(1)
    atype = config.get('structure', 'atype')
    atype = [a for a in atype.split()]    # list
    if not algo == 'EA-vc':
        natot = config.getint('structure', 'natot')
        if natot <= 0:
            logger.error('natot <= 0, check natot')
            raise SystemExit(1)
        nat = config.get('structure', 'nat')
        nat = [int(x) for x in nat.split()]    # character --> integer
        if not len(nat) == len(atype):
            logger.error('not len(nat) == len(atype), check atype and nat')
            raise SystemExit(1)
        if not sum(nat) == natot:
            logger.error('not sum(nat) == natot, check natot and nat')
            raise SystemExit(1)
    if algo =='EA-vc':
        ll_nat = config.get('structure', 'll_nat')
        ll_nat = [int(x) for x in ll_nat.split()]
        ul_nat = config.get('structure', 'ul_nat')
        ul_nat = [int(x) for x in ul_nat.split()]
        if not len(atype) == len(ll_nat) == len(ul_nat):
            logger.error('not len(atype) == len(ll_nat) == len(ul_nat), check ll_nat and ul_nat')
            raise SystemExit(1)
        for i in range(len(ll_nat)):
            if not 1 <= ll_nat[i] <= ul_nat[i]:
                logger.error(f'not 1 <= ll_nat[{i}] <= ul_nat[{i}], check ll_nat and ul_nat')
                raise SystemExit(1)
    # -- mol
    if struc_mode in ['mol', 'mol_bs']:
        mol_file = config.get('structure', 'mol_file')
        mol_file = [a for a in mol_file.split()]    # list
        nmol = config.get('structure', 'nmol')
        nmol = [int(x) for x in nmol.split()]    # character --> integer
        if not len(mol_file) == len(nmol):
            logger.error('not len(mol_file) == len(nmol)')
            raise SystemExit(1)
        try:
            timeout_mol = config.getfloat('structure', 'timeout_mol')
        except (configparser.NoOptionError, configparser.NoSectionError):
            timeout_mol = 120.0
        if timeout_mol <= 0:
            logger.error('timeout_mol must be positive')
            raise SystemExit(1)
        if struc_mode == 'mol_bs':
            # rot_mol
            try:
                rot_mol = config.get('structure', 'rot_mol')
            except (configparser.NoOptionError, configparser.NoSectionError):
                rot_mol = 'random_wyckoff'
            if rot_mol not in ['random', 'random_mol', 'random_wyckoff']:
                logger.error('rot_mol is wrong')
                raise SystemExit(1)
            # nrot
            try:
                nrot = config.getint('structure', 'nrot')
            except (configparser.NoOptionError, configparser.NoSectionError):
                nrot = 20
            if nrot <= 0:
                logger.error('nrot <=0, check nrot')
                raise SystemExit(1)
        else:
            rot_mol = None
            nrot = None
    else:
        mol_file = None
        nmol = None
        timeout_mol = 120.0
        rot_mol = None
        nrot = None
    # -- volume
    try:
        vol_mu = config.getfloat('structure', 'vol_mu')
    except (configparser.NoOptionError, configparser.NoSectionError):
        vol_mu = None
    if vol_mu is not None:
        if vol_mu <= 0.0:
            logger.error('vol_mu must be positive float')
            raise SystemExit(1)
    try:
        vol_sigma = config.getfloat('structure', 'vol_sigma')
    except (configparser.NoOptionError, configparser.NoSectionError):
        vol_sigma = None
    if vol_mu is not None:
        if vol_sigma is None:
            logger.error(f'check vol_mu: {vol_mu} and vol_sigma: {vol_sigma}')
            raise SystemExit(1)
    if vol_sigma is not None:
        if vol_sigma < 0.0:
            logger.error('vol_sigma must not be negative')
            raise SystemExit(1)
    try:
        vol_factor = config.get('structure', 'vol_factor')
        vol_factor = [float(x) for x in vol_factor.split()]    # char --> float
        if vol_factor[0] <= 0.0:
            logger.error('vol_factor must be positive')
            raise SystemExit(1)
        if len(vol_factor) == 1:
            vol_factor = vol_factor * 2    # [0.8] --> [0.8, 0.8]
        if len(vol_factor) == 2:
            if vol_factor[0] > vol_factor[1]:
                logger.error('check: vol_factor[0] < vol_factor[1]')
                raise SystemExit(1)
        else:
            logger.error('len(vol_factor) must be 1 or 2')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        vol_factor = [1.1, 1.1]
    try:
        maxcnt = config.getint('structure', 'maxcnt')
    except (configparser.NoOptionError, configparser.NoSectionError):
        maxcnt = 50
    if maxcnt < 0:
        logger.error('maxcnt must be positive int')
        raise SystemExit(1)
    try:
        symprec = config.getfloat('structure', 'symprec')
    except (configparser.NoOptionError, configparser.NoSectionError):
        symprec = 0.01
    if symprec < 0.0:
        logger.error('symprec must be positive float')
        raise SystemExit(1)
    try:
        spgnum = config.get('structure', 'spgnum')
    except (configparser.NoOptionError, configparser.NoSectionError):
        spgnum = 'all'
    if spgnum == '0':
        if struc_mode in ['mol', 'mol_bs']:
            logger.error('spgnum = 0 is not allow when struc_mode is mol or mol_bs')
            raise SystemExit(1)
        spgnum = 0
    elif spgnum == 'all':
        pass
    else:
        spgnum = _spglist(spgnum)
    try:
        use_find_wy = config.getboolean('structure', 'use_find_wy')
    except (configparser.NoOptionError, configparser.NoSectionError):
        use_find_wy = False
    if use_find_wy:
        if not struc_mode == 'crystal':
            logger.error('find_wy can be use if struc_mode is crystal')
            raise SystemExit(1)
    try:
        fwpath = config.get('structure', 'fwpath')
    except (configparser.NoOptionError, configparser.NoSectionError):
        fwpath = None
    # ------ mindist
    try:
        mindist = []
        for i in range(len(atype)):
            tmp = config.get('structure', f'mindist_{i+1}')
            tmp = [float(x) for x in tmp.split()]    # character --> float
            if not len(tmp) == len(atype):
                logger.error(f'not len(mindist_{i+1}) == len(atype)')
                raise SystemExit(1)
            mindist.append(tmp)
        # -- check symmetric matrix
        for i in range(len(mindist)):
            for j in range(len(mindist)):
                if i < j:
                    if not mindist[i][j] == mindist[j][i]:
                        logger.error(f'mindist is not symmetric. ({i}, {j}) -->'
                                         f' {mindist[i][j]}, ({j}, {i}) --> {mindist[j][i]}')
                        raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        mindist = None
    # ------ mindist_factor
    try:
        mindist_factor = config.getfloat('structure', 'mindist_factor')
        if mindist_factor <= 0.0:
            logger.error('mindist_factor must be positive')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        mindist_factor = 1.0
    # ------ mindist_mol_bs
    if struc_mode == 'mol_bs':
        try:
            mindist_mol_bs = []
            for i in range(len(mol_file)):
                tmp = config.get('structure', f'mindist_mol_bs_{i+1}')
                tmp = [float(x) for x in tmp.split()]    # character --> float
                if not len(tmp) == len(mol_file):
                    logger.error(f'not len(mindist_mol_bs_{i+1}) == len(mol_file)')
                    raise SystemExit(1)
                mindist_mol_bs.append(tmp)
            # -- check symmetric matrix
            for i in range(len(mindist_mol_bs)):
                for j in range(len(mindist_mol_bs)):
                    if i < j:
                        if not mindist_mol_bs[i][j] == mindist_mol_bs[j][i]:
                            logger.error(f'mindist_mol_bs is not symmetric. ({i}, {j}) -->'
                                             f' {mindist_mol_bs[i][j]}, ({j}, {i}) --> {mindist_mol_bs[j][i]}')
                            raise SystemExit(1)
        except (configparser.NoOptionError, configparser.NoSectionError):
            mindist_mol_bs = None
        # ------ mindist_mol_bs_factor
        try:
            mindist_mol_bs_factor = config.getfloat('structure', 'mindist_mol_bs_factor')
            if mindist_mol_bs_factor <= 0.0:
                logger.error('mindist_mol_bs_factor must be positive')
                raise SystemExit(1)
        except (configparser.NoOptionError, configparser.NoSectionError):
            mindist_mol_bs_factor = 1.0
    else:
        mindist_mol_bs = None
        mindist_mol_bs_factor = 1.0
    # ------ spgnum == 0 or use_find_wy
    minlen = None
    maxlen = None
    dangle = None
    if spgnum == 0 or use_find_wy:
        # -- read input variables
        minlen = config.getfloat('structure', 'minlen')
        maxlen = config.getfloat('structure', 'maxlen')
        dangle = config.getfloat('structure', 'dangle')
        if minlen <= 0.0:
            logger.error('minlen must be positive')
            raise SystemExit(1)
        if minlen > maxlen:
            logger.error('minlen > maxlen')
            raise SystemExit(1)
        if dangle <= 0.0:
            logger.error('dangle < 0.0, dangle must be positive')
            raise SystemExit(1)

    # ---------- option
    # ------ global declaration
    global stop_chkpt
    global load_struc_flag, stop_next_struc, recalc
    global append_struc_ea
    global energy_step_flag, struc_step_flag
    global force_step_flag, stress_step_flag

    # ------ read intput variables
    try:
        stop_chkpt = config.getint('option', 'stop_chkpt')
    except (configparser.NoOptionError, configparser.NoSectionError):
        stop_chkpt = 0
    try:
        load_struc_flag = config.getboolean('option', 'load_struc_flag')
    except (configparser.NoOptionError, configparser.NoSectionError):
        load_struc_flag = False
    try:
        stop_next_struc = config.getboolean('option', 'stop_next_struc')
    except (configparser.NoOptionError, configparser.NoSectionError):
        stop_next_struc = False
    try:
        recalc = config.get('option', 'recalc')
        recalc = [int(x) for x in recalc.split()]    # character --> integer
    except (configparser.NoOptionError, configparser.NoSectionError):
        recalc = []
    if recalc:
        for i in recalc:
            if not 0 <= i < tot_struc:
                logger.error('recalc must be non-negative int'
                                 ' and less than tot_struc')
                raise SystemExit(1)
    try:
        append_struc_ea = config.getboolean('option', 'append_struc_ea')
    except (configparser.NoOptionError, configparser.NoSectionError):
        append_struc_ea = False
    try:
        energy_step_flag = config.getboolean('option', 'energy_step_flag')
        if calc_code in ['LAMMPS', 'OMX', 'ASE', 'ext']:
            logger.error('energy_step_flag: only VASP, QE, and soiap for now')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        energy_step_flag = False
    try:
        struc_step_flag = config.getboolean('option', 'struc_step_flag')
        if calc_code in ['LAMMPS', 'OMX', 'ASE', 'ext']:
            logger.error('struc_step_flag: only VASP, QE, and soiap for now')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        struc_step_flag = False
    try:
        force_step_flag = config.getboolean('option', 'force_step_flag')
        if calc_code in ['LAMMPS', 'OMX', 'ASE', 'ext']:
            logger.error('force_step_flag: only VASP, QE, and soiap for now')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        force_step_flag = False
    if algo == 'LAQA':
        force_step_flag = True
    try:
        stress_step_flag = config.getboolean('option', 'stress_step_flag')
        if calc_code in ['LAMMPS', 'OMX', 'ASE', 'ext']:
            logger.error('stress_step_flag: only VASP, QE, and soiap for now')
            raise SystemExit(1)
    except (configparser.NoOptionError, configparser.NoSectionError):
        stress_step_flag = False
    if algo == 'LAQA':
        stress_step_flag = True

    # ---------- BO
    if algo == 'BO':
        # ------ global declaration
        global nselect_bo, score, num_rand_basis, cdev, dscrpt
        global fppath, fp_rmin, fp_rmax, fp_npoints, fp_sigma
        global max_select_bo, manual_select_bo, emax_bo, emin_bo
        # ------ read intput variables
        nselect_bo = config.getint('BO', 'nselect_bo')
        if nselect_bo <= 0:
            logger.error('nselect_bo <= 0, check nselect_bo')
            raise SystemExit(1)
        elif tot_struc < nselect_bo:
            logger.error('tot_struc < nselect_bo, check nselect_bo')
            raise SystemExit(1)
        score = config.get('BO', 'score')
        if score == 'TS' or score == 'EI' or score == 'PI':
            pass
        else:
            logger.error('score must be TS, EI, or PI, check score')
            raise SystemExit(1)
        try:
            num_rand_basis = config.getint('BO', 'num_rand_basis')
        except configparser.NoOptionError:
            num_rand_basis = 0
        try:
            cdev = config.getfloat('BO', 'cdev')
        except configparser.NoOptionError:
            cdev = 0.001
        dscrpt = config.get('BO', 'dscrpt')
        if dscrpt == 'FP':
            try:
                fppath = config.get('BO', 'fppath')
            except (configparser.NoOptionError, configparser.NoSectionError):
                fppath = None
            # -- check cal_fingerprint executable file
            fppath = utility.check_fppath(fppath)
        else:
            logger.error('Now FP only')
            raise SystemExit(1)
        # -- parameters for f-fingerprint
        try:
            fp_rmin = config.getfloat('BO', 'fp_rmin')
        except configparser.NoOptionError:
            fp_rmin = 0.5
        try:
            fp_rmax = config.getfloat('BO', 'fp_rmax')
        except configparser.NoOptionError:
            fp_rmax = 5.0
        if fp_rmin < 0.0:
            logger.error('fp_rmin < 0, check fp_rmin')
            raise SystemExit(1)
        if fp_rmax < fp_rmin:
            logger.error('fp_rmax < fp_rmin, check fp_rmin and fp_rmax')
            raise SystemExit(1)
        try:
            fp_npoints = config.getint('BO', 'fp_npoints')
        except configparser.NoOptionError:
            fp_npoints = 20
        if fp_npoints <= 0:
            logger.error('fp_npoints <= 0, check fp_npoints')
            raise SystemExit(1)
        try:
            fp_sigma = config.getfloat('BO', 'fp_sigma')
        except configparser.NoOptionError:
            fp_sigma = 1.0
        if fp_sigma < 0:
            logger.error('fp_sigma < 0, check fp_sigma')
            raise SystemExit(1)
        # -- BO option
        try:
            max_select_bo = config.getint('BO', 'max_select_bo')
        except configparser.NoOptionError:
            max_select_bo = 0
        if max_select_bo < 0:
            logger.error('max_select_bo must be non-negative int')
            raise SystemExit(1)
        try:
            manual_select_bo = config.get('BO', 'manual_select_bo')
            manual_select_bo = [int(x) for x in manual_select_bo.split()]
        except configparser.NoOptionError:
            manual_select_bo = []
        if manual_select_bo:
            for i in manual_select_bo:
                if not 0 <= i < tot_struc:
                    logger.error('manual_select_bo must be'
                                     ' non-negative int'
                                     ' and less than tot_struc')
                    raise SystemExit(1)
        try:
            emax_bo = config.getfloat('BO', 'emax_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            emax_bo = None
        try:
            emin_bo = config.getfloat('BO', 'emin_bo')
        except (configparser.NoOptionError, configparser.NoSectionError):
            emin_bo = None
        if emax_bo is not None and emin_bo is not None:
            if emin_bo > emax_bo:
                logger.error('emax_bo < emin_bo, check emax_bo and emin_bo')
                raise SystemExit(1)

    # ---------- LAQA
    if algo == 'LAQA':
        # ------ global declaration
        global nselect_laqa, wf, ws
        # ------ read intput variables
        nselect_laqa = config.getint('LAQA', 'nselect_laqa')
        try:
            wf = config.getfloat('LAQA', 'wf')
        except configparser.NoOptionError:
            wf = 0.1
        try:
            ws = config.getfloat('LAQA', 'ws')
        except configparser.NoOptionError:
            ws = 10.0

    # ---------- EA
    if algo in ['EA', 'EA-vc'] or append_struc_ea:
        # ------ global declaration
        global n_pop, n_crsov, n_perm, n_strain, n_rand, n_elite
        global fit_reverse, n_fittest
        global slct_func, t_size, a_rlt, b_rlt
        global crs_lat, nat_diff_tole, ntimes, sigma_st,  maxcnt_ea
        global maxgen_ea, emax_ea, emin_ea
        global end_point
        if struc_mode in ['mol', 'mol_bs']:
            global n_rotation, mindist_mol_ea, rot_max_angle, protect_mol_struc
        # global restart_gen
        # ------ read intput variables
        # -- number of structures
        n_pop = config.getint('EA', 'n_pop')
        if n_pop <= 0:
            logger.error('n_pop must be positive int')
            raise SystemExit(1)
        n_crsov = config.getint('EA', 'n_crsov')
        if n_crsov < 0:
            logger.error('n_crsov must be zero or positive int')
            raise SystemExit(1)
        n_perm = config.getint('EA', 'n_perm')
        if n_perm < 0:
            logger.error('n_perm must be zero or positive int')
            raise SystemExit(1)
        if n_perm != 0 and len(atype) == 1:
            logger.error('When the number of atom type is 1,'
                             ' n_perm must be 0')
            raise SystemExit(1)
        n_strain = config.getint('EA', 'n_strain')
        if n_strain < 0:
            logger.error('n_strain must be zero or positive int')
            raise SystemExit(1)
        n_rand = config.getint('EA', 'n_rand')
        if n_rand < 0:
            logger.error('n_rand must be zero or positive int')
            raise SystemExit(1)
        if struc_mode not in ['mol', 'mol_bs']:
            if n_crsov + n_perm + n_strain + n_rand != n_pop:
                logger.error('n_crsov + n_perm + n_strain + n_rand'
                                 ' must be n_pop')
                raise SystemExit(1)
        if struc_mode in ['mol', 'mol_bs']:
            n_rotation = config.getint('EA', 'n_rotation')
            if n_rotation < 0:
                logger.error('n_rotation must be zero or positive int')
                raise SystemExit(1)
            if n_crsov + n_perm + n_strain + n_rand + n_rotation != n_pop:
                logger.error('n_crsov + n_perm + n_strain + n_rand + n_rotation'
                                 ' must be n_pop')
                raise SystemExit(1)
        n_elite = config.getint('EA', 'n_elite')
        if n_elite < 0:
            logger.error('n_elite must be non-negative int')
            raise SystemExit(1)
        # -- n_fittest
        try:
            fit_reverse = config.getboolean('EA', 'fit_reverse')
        except configparser.NoOptionError:
            fit_reverse = False
        try:
            n_fittest = config.getint('EA', 'n_fittest')
        except configparser.NoOptionError:
            n_fittest = 0
        if n_fittest < 0:
            logger.error('n_fittest must be zero or positive int')
            raise SystemExit(1)
        # -- select function
        slct_func = config.get('EA', 'slct_func')
        if slct_func not in ['TNM', 'RLT']:
            logger.error('slct_func must be TNM or RLT')
            raise SystemExit(1)
        if slct_func == 'TNM':
            try:
                t_size = config.getint('EA', 't_size')
            except configparser.NoOptionError:
                t_size = 3
            if t_size < 2:
                logger.error('t_size must be greater than or equal to 2')
                raise SystemExit(1)
        elif slct_func == 'RLT':
            try:
                a_rlt = config.getfloat('EA', 'a_rlt')
            except configparser.NoOptionError:
                a_rlt = 10.0
            try:
                b_rlt = config.getfloat('EA', 'b_rlt')
            except configparser.NoOptionError:
                b_rlt = 1.0
            if not 0 < b_rlt < a_rlt:
                logger.error('must be 0 < b_rlt < a_rlt')
                raise SystemExit(1)
        # -- crossover
        try:
            crs_lat = config.get('EA', 'crs_lat')
        except configparser.NoOptionError:
            crs_lat = 'equal'
        if crs_lat not in ['equal', 'random']:
            logger.error('crs_lat must be equal or random')
            raise SystemExit(1)
        try:
            nat_diff_tole = config.getint('EA', 'nat_diff_tole')
        except configparser.NoOptionError:
            nat_diff_tole = 4
        if nat_diff_tole < 0:
            logger.error('nat_diff_tole must be nen-negative int')
            raise SystemExit(1)
        # -- permutation
        try:
            ntimes = config.getint('EA', 'ntimes')
        except configparser.NoOptionError:
            ntimes = 1
        if ntimes <= 0:
            logger.error('ntimes must be positive int')
            raise SystemExit(1)
        try:
            sigma_st = config.getfloat('EA', 'sigma_st')
        except configparser.NoOptionError:
            sigma_st = 0.5
        if sigma_st <= 0:
            logger.error('simga_st must be positive float')
            raise SystemExit(1)
        if struc_mode in ['mol', 'mol_bs']:
            # -- strain
            try:
                protect_mol_struc = config.getboolean('EA', 'protect_mol_struc')
            except configparser.NoOptionError:
                protect_mol_struc = False
            # -- rotation
            try:
                rot_max_angle = config.getint('EA', 'rot_max_angle')
            except configparser.NoOptionError:
                rot_max_angle = 360
            if rot_max_angle <= 0:
                logger.error('rot_max_angle must be positive int')
                raise SystemExit(1)
        # -- common
        if struc_mode in ['mol', 'mol_bs']:
            mindist_mol_ea = []
            for i in range(len(mol_file)):
                tmp = config.get('EA', f'mindist_mol_ea_{i+1}')
                tmp = [float(x) for x in tmp.split()]    # character --> float
                if not len(tmp) == len(mol_file):
                    logger.error(f'not len(mindist_mol_ea_{i+1}) == len(mol_file)')
                    raise SystemExit(1)
                mindist_mol_ea.append(tmp)
            # check symmetric matrix
            for i in range(len(mindist_mol_ea)):
                for j in range(len(mindist_mol_ea)):
                    if i < j:
                        if not mindist_mol_ea[i][j] == mindist_mol_ea[j][i]:
                            logger.error(f'mindist_mol is not symmetric. ({i}, {j}) -->'
                                             f' {mindist_mol_ea[i][j]}, ({j}, {i}) --> {mindist_mol_ea[j][i]}')
                            raise SystemExit(1)
        try:
            maxcnt_ea = config.getint('EA', 'maxcnt_ea')
        except configparser.NoOptionError:
            maxcnt_ea = 50
        # -- EA option
        try:
            maxgen_ea = config.getint('EA', 'maxgen_ea')
        except configparser.NoOptionError:
            maxgen_ea = 0
        if maxgen_ea < 0:
            logger.error('maxgen_ea must be non-negative int')
            raise SystemExit(1)
        # # -- restart option
        # try:
        #     restart_gen = config.getint('EA', 'restart_gen')
        # except configparser.NoOptionError:
        #     restart_gen = 0
        try:
            emax_ea = config.getfloat('EA', 'emax_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            emax_ea = None
        try:
            emin_ea = config.getfloat('EA', 'emin_ea')
        except (configparser.NoOptionError, configparser.NoSectionError):
            emin_ea = None
        if emax_ea is not None and emin_ea is not None:
            if emin_ea > emax_ea:
                logger.error('emax_ea < emin_ea, check emax_ea and emin_ea')
                raise SystemExit(1)
        # -- EA-vc
        if algo == 'EA-vc':
            end_point = config.get('EA','end_point')
            end_point = [float(x) for x in end_point.split()]
            if not len(end_point) == len(atype):
                logger.error('len(end_point) == len(atype), check end_point')
                raise SystemExit(1)
    # ---------- global declaration for comman part in calc_code
    global kppvol, kpt_flag, force_gamma

    # ---------- VASP
    if calc_code == 'VASP':
        # ------ read intput variables
        kpt_flag = True
        kppvol = config.get('VASP', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            logger.error('not len(kppvol) == nstage,'
                             ' check kppvol and nstage')
            raise SystemExit(1)
        try:
            force_gamma = config.getboolean('VASP', 'force_gamma')
        except configparser.NoOptionError:
            force_gamma = False

    # ---------- QE
    elif calc_code == 'QE':
        # ------ global declaration
        global qe_infile, qe_outfile, pv_term
        # ------ read intput variables
        kpt_flag = True
        qe_infile = config.get('QE', 'qe_infile')
        qe_outfile = config.get('QE', 'qe_outfile')
        kppvol = config.get('QE', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            logger.error('not len(kppvol) == nstage,'
                             ' check kppvol and nstage')
            raise SystemExit(1)
        try:
            force_gamma = config.getboolean('QE', 'force_gamma')
        except configparser.NoOptionError:
            force_gamma = False
        try:
            pv_term = config.getboolean('QE', 'pv_term')
        except configparser.NoOptionError:
            pv_term = False
        if energy_step_flag and pv_term:
            logger.error('cannot parse energy_step with pv_term yet. use pv_term = False')

    # ---------- OpenMX
    elif calc_code == 'OMX':
        # ------ global declaration
        global OMX_infile, OMX_outfile
        global upSpin, downSpin
        upSpin   = {}
        downSpin = {}
        # ------ read intput variables
        kpt_flag   = True
        OMX_infile  = config.get('OMX', 'OMX_infile')
        OMX_outfile = config.get('OMX', 'OMX_outfile')
        ValenceElec = config.get('OMX', 'ValenceElectrons')
        ValElecIn = ValenceElec.split()
        for i in range(0, len(ValElecIn), 3):
            upSpin[ValElecIn[i]]   = ValElecIn[i+1]
            downSpin[ValElecIn[i]] = ValElecIn[i+2]
        kppvol = config.get('OMX', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            logger.error('not len(kppvol) == nstage,'
                            ' check kppvol and nstage')
            raise SystemExit(1)
        try:
            force_gamma = config.getboolean('OMX', 'force_gamma')
        except configparser.NoOptionError:
            force_gamma = False

    # ---------- soiap
    elif calc_code == 'soiap':
        # ------ global declaration
        global soiap_infile, soiap_outfile, soiap_cif
        # ------ read intput variables
        soiap_infile = config.get('soiap', 'soiap_infile')
        soiap_outfile = config.get('soiap', 'soiap_outfile')
        soiap_cif = config.get('soiap', 'soiap_cif')
        kpt_flag = False
        force_gamma = False

    # ---------- lammps
    elif calc_code == 'LAMMPS':
        # ------ global declaration
        global lammps_infile, lammps_outfile, lammps_potential, lammps_data
        # ------ read intput variables
        lammps_infile = config.get('LAMMPS', 'lammps_infile')
        lammps_outfile = config.get('LAMMPS', 'lammps_outfile')
        try:
            lammps_potential = config.get('LAMMPS', 'lammps_potential')
            lammps_potential = lammps_potential.split()
        except configparser.NoOptionError:
            lammps_potential = None
        lammps_data = config.get('LAMMPS', 'lammps_data')
        kpt_flag = False
        force_gamma = False

    # ---------- ase
    elif calc_code == 'ASE':
        # ------ global declaration
        global ase_python
        # ------ read input variables
        ase_python = config.get('ASE', 'ase_python')
        kpt_flag = False
        force_gamma = False

    # ---------- ext
    elif calc_code == 'ext':
        kpt_flag = False


def _spglist(spgnum):
    tmpspg = []
    for c in spgnum.split():
        if '-' in c:
            if not len(c.split('-')) == 2:
                logger.error('Wrong input in spgnum. ')
                raise SystemExit(1)
            istart = int(c.split('-')[0])
            iend = int(c.split('-')[1])+1
            if istart < 0 or 230 < istart:
                logger.error('spgnum must be 1 -- 230')
                raise SystemExit(1)
            if iend < 0 or 231 < iend:
                logger.error('spgnum must be 1 -- 230')
                raise SystemExit(1)
            for i in range(istart, iend):
                if i not in tmpspg:
                    tmpspg.append(i)
        else:
            if int(c) < 0 or 230 < int(c):
                logger.error('spgnum must be 1 -- 230')
                raise SystemExit(1)
            if not int(c) in tmpspg:
                tmpspg += [int(c)]
    return tmpspg


def save_stat(stat):    # only 1st run
    logger.info('Save input data in cryspy.stat')
    # ---------- basic
    stat.set('basic', 'algo', '{}'.format(algo))
    stat.set('basic', 'calc_code', '{}'.format(calc_code))
    stat.set('basic', 'tot_struc', '{}'.format(tot_struc))
    if not calc_code == 'ext':
        stat.set('basic', 'nstage', '{}'.format(nstage))
        stat.set('basic', 'njob', '{}'.format(njob))
        stat.set('basic', 'jobcmd', '{}'.format(jobcmd))
        stat.set('basic', 'jobfile', '{}'.format(jobfile))

    # ---------- structure
    stat.set('structure', 'struc_mode', '{}'.format(struc_mode))
    stat.set('structure', 'atype', '{}'.format(' '.join(a for a in atype)))
    if not algo == 'EA-vc':
        stat.set('structure', 'natot', '{}'.format(natot))
        stat.set('structure', 'nat', '{}'.format(' '.join(str(b) for b in nat)))
    if algo == 'EA-vc':
        stat.set('structure', 'll_nat', '{}'.format(' '.join(str(b) for b in ll_nat)))
        stat.set('structure', 'ul_nat', '{}'.format(' '.join(str(b) for b in ul_nat)))
    if mol_file is None:
        stat.set('structure', 'mol_file', '{}'.format(mol_file))
    else:
        stat.set('structure', 'mol_file', '{}'.format(' '.join(a for a in mol_file)))
    if nmol is None:
        stat.set('structure', 'nmol', '{}'.format(nmol))
    else:
        stat.set('structure', 'nmol', '{}'.format(' '.join(str(b) for b in nmol)))
    stat.set('structure', 'timeout_mol', '{}'.format(timeout_mol))
    stat.set('structure', 'rot_mol', '{}'.format(rot_mol))
    stat.set('structure', 'nrot', '{}'.format(nrot))
    stat.set('structure', 'vol_factor', '{}'.format(' '.join(str(b) for b in vol_factor)))
    stat.set('structure', 'vol_mu', '{}'.format(vol_mu))
    stat.set('structure', 'vol_sigma', '{}'.format(vol_sigma))
    if mindist is None:
        stat.set('structure', 'mindist', '{}'.format(mindist))
    else:
        for i in range(len(atype)):
            stat.set('structure', 'mindist_{}'.format(i+1),
                     '{}'.format(' '.join(str(c) for c in mindist[i])))
    stat.set('structure', 'mindist_factor', '{}'.format(mindist_factor))
    if mindist_mol_bs is None:
        stat.set('structure', 'mindist_mol_bs', '{}'.format(mindist_mol_bs))
    else:
        for i in range(len(mol_file)):
            stat.set('structure', 'mindist_mol_bs_{}'.format(i+1),
                     '{}'.format(' '.join(str(c) for c in mindist_mol_bs[i])))
    stat.set('structure', 'mindist_mol_bs_factor', '{}'.format(mindist_mol_bs_factor))

    stat.set('structure', 'maxcnt', '{}'.format(maxcnt))
    stat.set('structure', 'symprec', '{}'.format(symprec))
    if spgnum == 0 or spgnum == 'all':
        stat.set('structure', 'spgnum', '{}'.format(spgnum))
    else:
        stat.set('structure', 'spgnum',
                 '{}'.format(' '.join(str(d) for d in spgnum)))
    stat.set('structure', 'use_find_wy', '{}'.format(use_find_wy))
    stat.set('structure', 'fwpath', '{}'.format(fwpath))
    stat.set('structure', 'minlen', '{}'.format(minlen))
    stat.set('structure', 'maxlen', '{}'.format(maxlen))
    stat.set('structure', 'dangle', '{}'.format(dangle))

    # ---------- BO
    if algo == 'BO':
        stat.set('BO', 'nselect_bo', '{}'.format(nselect_bo))
        stat.set('BO', 'score', '{}'.format(score))
        stat.set('BO', 'num_rand_basis', '{}'.format(num_rand_basis))
        stat.set('BO', 'cdev', '{}'.format(cdev))
        stat.set('BO', 'dscrpt', '{}'.format(dscrpt))
        stat.set('BO', 'fppath', '{}'.format(fppath))
        stat.set('BO', 'fp_rmin', '{}'.format(fp_rmin))
        stat.set('BO', 'fp_rmax', '{}'.format(fp_rmax))
        stat.set('BO', 'fp_npoints', '{}'.format(fp_npoints))
        stat.set('BO', 'fp_sigma', '{}'.format(fp_sigma))
        stat.set('BO', 'max_select_bo', '{}'.format(max_select_bo))
        stat.set('BO', 'manual_select_bo', '{}'.format(
            ' '.join(str(x) for x in manual_select_bo)))
        stat.set('BO', 'emax_bo', '{}'.format(emax_bo))
        stat.set('BO', 'emin_bo', '{}'.format(emin_bo))

    # ---------- LAQA
    if algo == 'LAQA':
        stat.set('LAQA', 'nselect_laqa', '{}'.format(nselect_laqa))
        stat.set('LAQA', 'wf', '{}'.format(wf))
        stat.set('LAQA', 'ws', '{}'.format(ws))

    # ---------- EA
    elif algo in ['EA', 'EA-vc'] or append_struc_ea:
        stat.set('EA', 'n_pop', '{}'.format(n_pop))
        stat.set('EA', 'n_crsov', '{}'.format(n_crsov))
        stat.set('EA', 'n_perm', '{}'.format(n_perm))
        stat.set('EA', 'n_strain', '{}'.format(n_strain))
        stat.set('EA', 'n_rand', '{}'.format(n_rand))
        if struc_mode in ['mol', 'mol_bs']:
            stat.set('EA', 'n_rotation', '{}'.format(n_rotation))
        stat.set('EA', 'n_elite', '{}'.format(n_elite))
        stat.set('EA', 'fit_reverse', '{}'.format(fit_reverse))
        stat.set('EA', 'n_fittest', '{}'.format(n_fittest))
        stat.set('EA', 'slct_func', '{}'.format(slct_func))
        if slct_func == 'TNM':
            stat.set('EA', 't_size', '{}'.format(t_size))
        elif slct_func == 'RLT':
            stat.set('EA', 'a_rlt', '{}'.format(a_rlt))
            stat.set('EA', 'b_rlt', '{}'.format(b_rlt))
        stat.set('EA', 'crs_lat', '{}'.format(crs_lat))
        stat.set('EA', 'nat_diff_tole', '{}'.format(nat_diff_tole))
        stat.set('EA', 'ntimes', '{}'.format(ntimes))
        stat.set('EA', 'sigma_st', '{}'.format(sigma_st))
        if struc_mode in ['mol', 'mol_bs']:
            stat.set('EA', 'protect_mol_struc', '{}'.format(protect_mol_struc))
            stat.set('EA', 'rot_max_angle', '{}'.format(rot_max_angle))
            for i, x in enumerate(mindist_mol_ea):
                stat.set('EA', 'mindist_mol_ea_{}'.format(i+1),
                         '{}'.format(' '.join(str(c) for c in x)))
        stat.set('EA', 'maxcnt_ea', '{}'.format(maxcnt_ea))
        stat.set('EA', 'maxgen_ea', '{}'.format(maxgen_ea))
        stat.set('EA', 'emax_ea', '{}'.format(emax_ea))
        stat.set('EA', 'emin_ea', '{}'.format(emin_ea))
        if algo == 'EA-vc':
            stat.set('EA', 'end_point', '{}'.format(' '.join(str(b) for b in end_point)))

    # ---------- VASP
    if calc_code == 'VASP':
        stat.set('VASP', 'kppvol',
                 '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('VASP', 'force_gamma', '{}'.format(force_gamma))

    # ---------- QE
    if calc_code == 'QE':
        stat.set('QE', 'qe_infile', '{}'.format(qe_infile))
        stat.set('QE', 'qe_outfile', '{}'.format(qe_outfile))
        stat.set('QE', 'kppvol',
                 '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('QE', 'force_gamma', '{}'.format(force_gamma))
        stat.set('QE', 'pv_term', '{}'.format(pv_term))

    # ---------- OMX
    if calc_code == 'OMX':
        stat.set('OMX', 'OMX_infile', '{}'.format(OMX_infile))
        stat.set('OMX', 'OMX_outfile', '{}'.format(OMX_outfile))
        stat.set('OMX', 'kppvol',
                 '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('OMX', 'force_gamma', '{}'.format(force_gamma))

    # ---------- soiap
    if calc_code == 'soiap':
        stat.set('soiap', 'soiap_infile', '{}'.format(soiap_infile))
        stat.set('soiap', 'soiap_outfile', '{}'.format(soiap_outfile))
        stat.set('soiap', 'soiap_cif', '{}'.format(soiap_cif))

    # ---------- lammps
    if calc_code == 'LAMMPS':
        stat.set('LAMMPS', 'lammps_infile', '{}'.format(lammps_infile))
        stat.set('LAMMPS', 'lammps_outfile', '{}'.format(lammps_outfile))
        if lammps_potential is None:
            stat.set('LAMMPS', 'lammps_potential',
                     '{}'.format(lammps_potential))
        else:
            stat.set('LAMMPS', 'lammps_potential',
                     '{}'.format(' '.join(lammps_potential)))
        stat.set('LAMMPS', 'lammps_data', '{}'.format(lammps_data))

    # ---------- ASE
    if calc_code == 'ASE':
        stat.set('ASE', 'ase_python', '{}'.format(ase_python))

    # ---------- option
    stat.set('option', 'stop_chkpt', '{}'.format(stop_chkpt))
    stat.set('option', 'load_struc_flag', '{}'.format(load_struc_flag))
    stat.set('option', 'stop_next_struc', '{}'.format(stop_next_struc))
    stat.set('option', 'recalc', '{}'.format(' '.join(str(x) for x in recalc)))
    stat.set('option', 'append_struc_ea', '{}'.format(append_struc_ea))
    stat.set('option', 'energy_step_flag', '{}'.format(energy_step_flag))
    stat.set('option', 'struc_step_flag', '{}'.format(struc_step_flag))
    stat.set('option', 'force_step_flag', '{}'.format(force_step_flag))
    stat.set('option', 'stress_step_flag', '{}'.format(stress_step_flag))

    # ---------- write stat
    io_stat.write_stat(stat)


def diffinstat(stat):
    #
    # Do not use SystemExit in diffinstat!
    #
    logic_change = False

    # ---------- old input
    # ------ basic
    old_algo = stat.get('basic', 'algo')
    old_calc_code = stat.get('basic', 'calc_code')
    old_tot_struc = stat.getint('basic', 'tot_struc')
    if not calc_code == 'ext':
        old_nstage = stat.getint('basic', 'nstage')
        old_njob = stat.getint('basic', 'njob')
        old_jobcmd = stat.get('basic', 'jobcmd')
        old_jobfile = stat.get('basic', 'jobfile')

    # ------ structure
    old_struc_mode = stat.get('structure', 'struc_mode')
    old_atype = stat.get('structure', 'atype')
    old_atype = [a for a in old_atype.split()]    # list
    if not algo == 'EA-vc':
        old_natot = stat.getint('structure', 'natot')
        old_nat = stat.get('structure', 'nat')
        old_nat = [int(x) for x in old_nat.split()]    # str --> int list
    if algo == 'EA-vc':
        old_ll_nat = stat.get('structure','ll_nat')
        old_ll_nat = [int(x) for x in old_ll_nat.split()]
        old_ul_nat = stat.get('structure','ul_nat')
        old_ul_nat = [int(x) for x in old_ul_nat.split()]
    old_mol_file = stat.get('structure', 'mol_file')
    if old_mol_file == 'None':
        old_mol_file = None    # character --> None
    else:
        old_mol_file = [a for a in old_mol_file.split()]    # list
    old_nmol = stat.get('structure', 'nmol')
    if old_nmol == 'None':
        old_nmol = None    # character --> None
    else:
        old_nmol = [int(x) for x in old_nmol.split()]    # str --> int list
    old_timeout_mol = stat.getfloat('structure', 'timeout_mol')
    old_rot_mol = stat.get('structure', 'rot_mol')
    if old_rot_mol == 'None':
        old_rot_mol = None    # character --> None
    old_nrot = stat.get('structure', 'nrot')
    if old_nrot == 'None':
        old_nrot = None    # character --> None
    else:
        old_nrot = int(old_nrot)    # character --> int
    old_vol_factor = stat.get('structure', 'vol_factor')
    old_vol_factor = [float(x) for x in old_vol_factor.split()]    # str --> float list
    old_vol_mu = stat.get('structure', 'vol_mu')
    if old_vol_mu == 'None':
        old_vol_mu = None    # character --> None
    else:
        old_vol_mu = float(old_vol_mu)    # character --> float
    old_vol_sigma = stat.get('structure', 'vol_sigma')
    if old_vol_sigma == 'None':
        old_vol_sigma = None    # character --> None
    else:
        old_vol_sigma = float(old_vol_sigma)    # character --> float
    try:    # case: None
        old_mindist = stat.get('structure', 'mindist')
        if old_mindist == 'None':
            old_mindist = None    # character --> None
    except (configparser.NoOptionError, configparser.NoSectionError):
        old_mindist = []
        for i in range(len(atype)):
            tmp = stat.get('structure', f'mindist_{i+1}')
            tmp = [float(x) for x in tmp.split()]    # character --> float
            old_mindist.append(tmp)
    old_mindist_factor = stat.getfloat('structure', 'mindist_factor')
    try:    # case: None
        old_mindist_mol_bs = stat.get('structure', 'mindist_mol_bs')
        if old_mindist_mol_bs == 'None':
            old_mindist_mol_bs = None    # character --> None
    except (configparser.NoOptionError, configparser.NoSectionError):
        old_mindist_mol_bs = []
        for i in range(len(mol_file)):
            tmp = stat.get('structure', f'mindist_mol_bs_{i+1}')
            tmp = [float(x) for x in tmp.split()]    # character --> float
            old_mindist_mol_bs.append(tmp)
    old_mindist_mol_bs_factor = stat.getfloat('structure', 'mindist_mol_bs_factor')
    old_maxcnt = stat.getint('structure', 'maxcnt')
    old_symprec = stat.getfloat('structure', 'symprec')
    old_spgnum = stat.get('structure', 'spgnum')
    if old_spgnum == '0':
        old_spgnum = 0
    elif not old_spgnum == 'all':
        old_spgnum = [int(x) for x in old_spgnum.split()]    # int list
    old_use_find_wy = stat.getboolean('structure', 'use_find_wy')
    old_fwpath = stat.get('structure', 'fwpath')
    if old_fwpath == 'None':
        old_fwpath = None    # character --> None
    old_minlen = stat.get('structure', 'minlen')
    if old_minlen == 'None':
        old_minlen = None    # character --> None
    else:
        old_minlen = float(old_minlen)    # character --> float
    old_maxlen = stat.get('structure', 'maxlen')
    if old_maxlen == 'None':
        old_maxlen = None    # character --> None
    else:
        old_maxlen = float(old_maxlen)    # character --> float
    old_dangle = stat.get('structure', 'dangle')
    if old_dangle == 'None':
        old_dangle = None    # character --> None
    else:
        old_dangle = float(old_dangle)    # character --> float

    # ------ BO
    if old_algo == 'BO':
        old_nselect_bo = stat.getint('BO', 'nselect_bo')
        old_score = stat.get('BO', 'score')
        old_num_rand_basis = stat.getint('BO', 'num_rand_basis')
        old_cdev = stat.getfloat('BO', 'cdev')
        old_dscrpt = stat.get('BO', 'dscrpt')
        old_fppath = stat.get('BO', 'fppath')
        if old_fppath == 'None':
            old_fppath = None    # character --> None
        old_fp_rmin = stat.getfloat('BO', 'fp_rmin')
        old_fp_rmax = stat.getfloat('BO', 'fp_rmax')
        old_fp_npoints = stat.getint('BO', 'fp_npoints')
        old_fp_sigma = stat.getfloat('BO', 'fp_sigma')
        old_max_select_bo = stat.getint('BO', 'max_select_bo')
        old_manual_select_bo = stat.get('BO', 'manual_select_bo')
        old_manual_select_bo = [int(x) for x in old_manual_select_bo.split()]
        old_emax_bo = stat.get('BO', 'emax_bo')
        if old_emax_bo == 'None':
            old_emax_bo = None    # char --> None
        else:
            old_emax_bo = float(old_emax_bo)    # char --> float
        old_emin_bo = stat.get('BO', 'emin_bo')
        if old_emin_bo == 'None':
            old_emin_bo = None    # char --> None
        else:
            old_emin_bo = float(old_emin_bo)    # char --> float

    # ------ LAQA
    if old_algo == 'LAQA':
        old_nselect_laqa = stat.getint('LAQA', 'nselect_laqa')
        old_wf = stat.getfloat('LAQA', 'wf')
        old_ws = stat.getfloat('LAQA', 'ws')


    # ------ EA
    if old_algo in ['EA', 'EA-vc']:
        old_n_pop = stat.getint('EA', 'n_pop')
        old_n_crsov = stat.getint('EA', 'n_crsov')
        old_n_perm = stat.getint('EA', 'n_perm')
        old_n_strain = stat.getint('EA', 'n_strain')
        old_n_rand = stat.getint('EA', 'n_rand')
        if struc_mode in ['mol', 'mol_bs']:
            old_n_rotation = stat.getint('EA', 'n_rotation')
        old_n_elite = stat.getint('EA', 'n_elite')
        old_fit_reverse = stat.getboolean('EA', 'fit_reverse')
        old_n_fittest = stat.getint('EA', 'n_fittest')
        old_slct_func = stat.get('EA', 'slct_func')
        if old_slct_func == 'TNM':
            old_t_size = stat.getint('EA', 't_size')
        elif old_slct_func == 'RLT':
            old_a_rlt = stat.getfloat('EA', 'a_rlt')
            old_b_rlt = stat.getfloat('EA', 'b_rlt')
        old_crs_lat = stat.get('EA', 'crs_lat')
        old_nat_diff_tole = stat.getint('EA', 'nat_diff_tole')
        old_ntimes = stat.getint('EA', 'ntimes')
        old_sigma_st = stat.getfloat('EA', 'sigma_st')
        if struc_mode in ['mol', 'mol_bs']:
            old_protect_mol_struc = stat.getboolean('EA', 'protect_mol_struc')
            old_rot_max_angle = stat.getint('EA', 'rot_max_angle')
            old_mindist_mol_ea = []
            for i in range(len(mol_file)):
                tmp = stat.get('EA', f'mindist_mol_ea_{i+1}')
                tmp = [float(x) for x in tmp.split()]    # character --> float
                old_mindist_mol_ea.append(tmp)
        old_maxcnt_ea = stat.getint('EA', 'maxcnt_ea')
        old_maxgen_ea = stat.getint('EA', 'maxgen_ea')
        # old_restart_gen = stat.get('EA', 'restart_gen')
        old_emax_ea = stat.get('EA', 'emax_ea')
        if old_emax_ea == 'None':
            old_emax_ea = None    # char --> None
        else:
            old_emax_ea = float(old_emax_ea)    # char --> float
        old_emin_ea = stat.get('EA', 'emin_ea')
        if old_emin_ea == 'None':
            old_emin_ea = None    # char --> None
        else:
            old_emin_ea = float(old_emin_ea)    # char --> float
        if algo == 'EA-vc':
            old_end_point = stat.get('EA','end_point')
            old_end_point = [float(x) for x in old_end_point.split()]

    # ------ VASP
    if old_calc_code == 'VASP':
        old_kppvol = stat.get('VASP', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]    # int list
        old_force_gamma = stat.getboolean('VASP', 'force_gamma')

    # ------ QE
    if old_calc_code == 'QE':
        old_qe_infile = stat.get('QE', 'qe_infile')
        old_qe_outfile = stat.get('QE', 'qe_outfile')
        old_kppvol = stat.get('QE', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]  # int list
        old_force_gamma = stat.getboolean('QE', 'force_gamma')
        old_pv_term = stat.getboolean('QE', 'pv_term')

    if old_calc_code == 'OMX':
        old_OMX_infile = stat.get('OMX', 'OMX_infile')
        old_OMX_outfile = stat.get('OMX', 'OMX_outfile')
        old_kppvol = stat.get('OMX', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]  # int list
        old_force_gamma = stat.getboolean('OMX', 'force_gamma')

    # ------ soiap
    if old_calc_code == 'soiap':
        old_soiap_infile = stat.get('soiap', 'soiap_infile')
        old_soiap_outfile = stat.get('soiap', 'soiap_outfile')
        old_soiap_cif = stat.get('soiap', 'soiap_cif')

    # ------ lammps
    if old_calc_code == 'LAMMPS':
        old_lammps_infile = stat.get('LAMMPS', 'lammps_infile')
        old_lammps_outfile = stat.get('LAMMPS', 'lammps_outfile')
        old_lammps_potential = stat.get('LAMMPS', 'lammps_potential')
        if old_lammps_potential == 'None':    # 'None' is just character here
            old_lammps_potential = None
        else:
            old_lammps_potential = old_lammps_potential.split()    # str --> list
        old_lammps_data = stat.get('LAMMPS', 'lammps_data')

    # ------ ASE
    if old_calc_code == 'ASE':
        old_ase_python = stat.get('ASE', 'ase_python')

    # ------ option
    old_stop_chkpt = stat.getint('option', 'stop_chkpt')
    old_load_struc_flag = stat.getboolean('option', 'load_struc_flag')
    old_stop_next_struc = stat.getboolean('option', 'stop_next_struc')
    old_recalc = stat.get('option', 'recalc')
    old_recalc = [int(x) for x in old_recalc.split()]    # int list
    old_append_struc_ea = stat.getboolean('option', 'append_struc_ea')
    old_energy_step_flag = stat.getboolean('option', 'energy_step_flag')
    old_struc_step_flag = stat.getboolean('option', 'struc_step_flag')
    old_force_step_flag = stat.getboolean('option', 'force_step_flag')
    old_stress_step_flag = stat.getboolean('option', 'stress_step_flag')

    # ---------- check difference
    # ------ basic
    sec = 'basic'
    if not old_algo == algo:
        raise ValueError('Do not change algo')
    if not old_calc_code == calc_code:
        raise ValueError('Do not change calc code')
    if not old_tot_struc == tot_struc:
        if algo in ['EA', 'EA-vc']:
            raise ValueError('Do not change tot_struc in EA')
        diff_out('tot_struc', old_tot_struc, tot_struc)
        io_stat.set_input_common(stat, sec, 'tot_struc', tot_struc)
        logic_change = True
    if not calc_code == 'ext':
        if not old_nstage == nstage:
            diff_out('nstage', old_nstage, nstage)
            io_stat.set_input_common(stat, sec, 'nstage', nstage)
            logic_change = True
        if not old_njob == njob:
            diff_out('njob', old_njob, njob)
            io_stat.set_input_common(stat, sec, 'njob', njob)
            logic_change = True
        if not old_jobcmd == jobcmd:
            diff_out('jobcmd', old_jobcmd, jobcmd)
            io_stat.set_input_common(stat, sec, 'jobcmd', jobcmd)
            logic_change = True
        if not old_jobfile == jobfile:
            diff_out('jobfile', old_jobfile, jobfile)
            io_stat.set_input_common(stat, sec, 'jobfile', jobfile)
            logic_change = True

    # ------ structure
    sec = 'structure'
    if not old_struc_mode == struc_mode:
        if old_struc_mode in ['crystal', 'mol', 'mol_bs'] and struc_mode in ['crystal', 'mol', 'mol_bs']:
            diff_out('struc_mode', old_struc_mode, struc_mode)
            io_stat.set_input_common(stat, sec, 'struc_mode', struc_mode)
            logic_change = True
        else:
            raise ValueError('Do not change struc_mode: host')
    if not old_atype == atype:
        raise ValueError('Do not change atype')
    if not algo == 'EA-vc':
        if not old_natot == natot:
            raise ValueError('Do not change natot')
        if not old_nat == nat:
            raise ValueError('Do not change nat')
    if algo == 'EA-vc':
        if not old_ll_nat == ll_nat:
            diff_out('ll_nat', old_ll_nat, ll_nat)
            io_stat.set_input_common(stat, sec, 'll_nat', ll_nat)
            logic_change = True
        if not old_ul_nat == ul_nat:
            diff_out('ul_nat', old_ul_nat, ul_nat)
            io_stat.set_input_common(stat, sec, 'ul_nat', ul_nat)
            logic_change = True
    if not old_mol_file == mol_file:
        diff_out('mol_file', old_mol_file, mol_file)
        io_stat.set_input_common(stat, sec, 'mol_file', mol_file)
        logic_change = True
    if not old_nmol == nmol:
        if old_nmol is None or nmol is None:
            diff_out('nmol', old_nmol, nmol)
            io_stat.set_input_common(stat, sec, 'nmol', nmol)
            logic_change = True
        else:
            raise ValueError('Do not change nmol except for None')
    if not old_timeout_mol == timeout_mol:
        diff_out('timeout_mol', old_timeout_mol, timeout_mol)
        io_stat.set_input_common(stat, sec, 'timeout_mol', timeout_mol)
        logic_change = True
    if not old_rot_mol == rot_mol:
        diff_out('rot_mol', old_rot_mol, rot_mol)
        io_stat.set_input_common(stat, sec, 'rot_mol', rot_mol)
        logic_change = True
    if not old_nrot == nrot:
        diff_out('nrot', old_nrot, nrot)
        io_stat.set_input_common(stat, sec, 'nrot', nrot)
        logic_change = True
    if not old_vol_factor == vol_factor:
        diff_out('vol_factor', old_vol_factor, vol_factor)
        io_stat.set_input_common(stat, sec, 'vol_factor', '{}'.format(
                ' '.join(str(x) for x in vol_factor)))
        logic_change = True
    if not old_vol_mu == vol_mu:
        diff_out('vol_mu', old_vol_mu, vol_mu)
        io_stat.set_input_common(stat, sec, 'vol_mu', vol_mu)
        logic_change = True
    if not old_vol_sigma == vol_sigma:
        diff_out('vol_sigma', old_vol_sigma, vol_sigma)
        io_stat.set_input_common(stat, sec, 'vol_sigma', vol_sigma)
        logic_change = True
    if not old_mindist == mindist:
        diff_out('mindist', old_mindist, mindist)
        # -- case: old_mindist = None, mindist = []
        if old_mindist is None:
            stat.remove_option('structure', 'mindist')    # clear mindist
            for i in range(len(atype)):               # add mindist_?
                io_stat.set_input_common(stat, sec, 'mindist_{}'.format(i+1),
                                         '{}'.format(' '.join(
                                             str(x) for x in mindist[i])))
        # -- case: old_mindist = [], mindist = None
        elif mindist is None:
            for i in range(len(atype)):    # clear mindist_?
                stat.remove_option('structure', 'mindist_{}'.format(i+1))
            io_stat.set_input_common(stat, sec, 'mindist', mindist)    # add mindist
        # -- case: old_mindist = [], mindist = [], update list
        else:
            for i in range(len(atype)):
                io_stat.set_input_common(stat, sec, 'mindist_{}'.format(i+1),
                                         '{}'.format(' '.join(
                                             str(x) for x in mindist[i])))
        logic_change = True
    if not old_mindist_factor == mindist_factor:
        diff_out('mindist_factor', old_mindist_factor, mindist_factor)
        io_stat.set_input_common(stat, sec, 'mindist_factor', mindist_factor)
        logic_change = True
    if not old_mindist_mol_bs == mindist_mol_bs:
        diff_out('mindist_mol_bs', old_mindist_mol_bs, mindist_mol_bs)
        # -- case: old_mindist_mol_bs = None, mindist_mol_bs = []
        if old_mindist_mol_bs is None:
            stat.remove_option('structure', 'mindist_mol_bs')    # clear mindist_mol_bs
            for i in range(len(mol_file)):               # add mindist_mol_bs_?
                io_stat.set_input_common(stat, sec, 'mindist_mol_bs_{}'.format(i+1),
                                         '{}'.format(' '.join(
                                             str(x) for x in mindist_mol_bs[i])))
        # -- case: old_mindist_mol_bs = [], mindist_mol_bs = None
        elif mindist_mol_bs is None:
            for i in range(len(mol_file)):    # clear mindist_mol_bs_?
                stat.remove_option('structure', 'mindist_mol_bs_{}'.format(i+1))
            io_stat.set_input_common(stat, sec, 'mindist_mol_bs', mindist_mol_bs)    # add mindist_mol_bs
        # -- case: old_mindist_mol_bs = [], mindist_mol_bs = [], update list
        else:
            for i in range(len(mol_file)):
                io_stat.set_input_common(stat, sec, 'mindist_mol_bs_{}'.format(i+1),
                                         '{}'.format(' '.join(
                                             str(x) for x in mindist_mol_bs[i])))
        logic_change = True
    if not old_mindist_mol_bs_factor == mindist_mol_bs_factor:
        diff_out('mindist_mol_bs_factor', old_mindist_mol_bs_factor, mindist_mol_bs_factor)
        io_stat.set_input_common(stat, sec, 'mindist_mol_bs_factor', mindist_mol_bs_factor)
        logic_change = True
    if not old_maxcnt == maxcnt:
        diff_out('maxcnt', old_maxcnt, maxcnt)
        io_stat.set_input_common(stat, sec, 'maxcnt', maxcnt)
        logic_change = True
    if not old_symprec == symprec:
        diff_out('symprec', old_symprec, symprec)
        io_stat.set_input_common(stat, sec, 'symprec', symprec)
        logic_change = True
    if not old_spgnum == spgnum:
        diff_out('spgnum', old_spgnum, spgnum)
        if spgnum == 0 or spgnum == 'all':
            io_stat.set_input_common(stat, sec, 'spgnum', spgnum)
        else:
            io_stat.set_input_common(stat, sec, 'spgnum', '{}'.format(
                ' '.join(str(x) for x in spgnum)))
        logic_change = True
    if not old_use_find_wy == use_find_wy:
        diff_out('use_find_wy', old_use_find_wy, use_find_wy)
        io_stat.set_input_common(stat, sec, 'use_find_wy', use_find_wy)
        logic_change = True
    if not old_fwpath == fwpath:
        diff_out('fwpath', old_fwpath, fwpath)
        io_stat.set_input_common(stat, sec, 'fwpath', fwpath)
        logic_change = True
    if not old_minlen == minlen:
        diff_out('minlen', old_minlen, minlen)
        io_stat.set_input_common(stat, sec, 'minlen', minlen)
        logic_change = True
    if not old_maxlen == maxlen:
        diff_out('maxlen', old_maxlen, maxlen)
        io_stat.set_input_common(stat, sec, 'maxlen', maxlen)
        logic_change = True
    if not old_dangle == dangle:
        diff_out('dangle', old_dangle, dangle)
        io_stat.set_input_common(stat, sec, 'dangle', dangle)
        logic_change = True

    # ------ BO
    sec = 'BO'
    if algo == 'BO':
        if not old_nselect_bo == nselect_bo:
            diff_out('nselect_bo', old_nselect_bo, nselect_bo)
            io_stat.set_input_common(stat, sec, 'nselect_bo', nselect_bo)
            logic_change = True
        if not old_score == score:
            diff_out('score', old_score, score)
            io_stat.set_input_common(stat, sec, 'score', score)
            logic_change = True
        if not old_num_rand_basis == num_rand_basis:
            diff_out('num_rand_basis', old_num_rand_basis, num_rand_basis)
            io_stat.set_input_common(stat, sec, 'num_rand_basis', num_rand_basis)
            logic_change = True
        if not old_cdev == cdev:
            diff_out('cdev', old_cdev, cdev)
            io_stat.set_input_common(stat, sec, 'cdev', cdev)
            logic_change = True
        if not old_dscrpt == dscrpt:
            raise ValueError('Do not change dscrpt')
        if not old_fppath == fppath:
            diff_out('fppath', old_fppath, fppath)
            io_stat.set_input_common(stat, sec, 'fppath', fppath)
            logic_change = True
        if not old_fp_rmin == fp_rmin:
            raise ValueError('Do not change fp_rmin')
        if not old_fp_rmax == fp_rmax:
            raise ValueError('Do not change fp_rmax')
        if not old_fp_npoints == fp_npoints:
            raise ValueError('Do not change fp_npoints')
        if not old_fp_sigma == fp_sigma:
            raise ValueError('Do not change fp_sigma')
        if not old_max_select_bo == max_select_bo:
            diff_out('max_select_bo', old_max_select_bo, max_select_bo)
            io_stat.set_input_common(stat, sec, 'max_select_bo', max_select_bo)
            logic_change = True
        if not old_manual_select_bo == manual_select_bo:
            diff_out('manual_select_bo', old_manual_select_bo,
                     manual_select_bo)
            io_stat.set_input_common(stat, sec, 'manual_select_bo',
                                     '{}'.format(' '.join(
                                         str(x) for x in manual_select_bo)))
            logic_change = True
        if not old_emax_bo == emax_bo:
            diff_out('emax_bo', old_emax_bo, emax_bo)
            io_stat.set_input_common(stat, sec, 'emax_bo', emax_bo)
            logic_change = True
        if not old_emin_bo == emin_bo:
            diff_out('emin_bo', old_emin_bo, emin_bo)
            io_stat.set_input_common(stat, sec, 'emin_bo', emin_bo)
            logic_change = True

    # ------ LAQA
    sec = 'LAQA'
    if algo == 'LAQA':
        if not old_nselect_laqa == nselect_laqa:
            diff_out('nselect_laqa', old_nselect_laqa, nselect_laqa)
            io_stat.set_input_common(stat, sec, 'nselect_laqa', nselect_laqa)
            logic_change = True
        if not old_wf == wf:
            diff_out('wf', old_wf, wf)
            io_stat.set_input_common(stat, sec, 'wf', wf)
            logic_change = True
        if not old_ws == ws:
            diff_out('ws', old_ws, ws)
            io_stat.set_input_common(stat, sec, 'ws', ws)
            logic_change = True

    # ------ EA
    sec = 'EA'
    if algo in ['EA', 'EA-vc']:
        if not old_n_pop == n_pop:
            diff_out('n_pop', old_n_pop, n_pop)
            io_stat.set_input_common(stat, sec, 'n_pop', n_pop)
            logic_change = True
        if not old_n_crsov == n_crsov:
            diff_out('n_crsov', old_n_crsov, n_crsov)
            io_stat.set_input_common(stat, sec, 'n_crsov', n_crsov)
            logic_change = True
        if not old_n_perm == n_perm:
            diff_out('n_perm', old_n_perm, n_perm)
            io_stat.set_input_common(stat, sec, 'n_perm', n_perm)
            logic_change = True
        if not old_n_strain == n_strain:
            diff_out('n_strain', old_n_strain, n_strain)
            io_stat.set_input_common(stat, sec, 'n_strain', n_strain)
            logic_change = True
        if not old_n_rand == n_rand:
            diff_out('n_rand', old_n_rand, n_rand)
            io_stat.set_input_common(stat, sec, 'n_rand', n_rand)
            logic_change = True
        if struc_mode in ['mol', 'mol_bs']:
            if not old_n_rotation == n_rotation:
                diff_out('n_rotation', old_n_rotation, n_rotation)
                io_stat.set_input_common(stat, sec, 'n_rotation', n_rotation)
                logic_change = True
        if not old_n_elite == n_elite:
            diff_out('n_elite', old_n_elite, n_elite)
            io_stat.set_input_common(stat, sec, 'n_elite', n_elite)
            logic_change = True
        if not old_fit_reverse == fit_reverse:
            raise ValueError('Do not change fit_reverse')
        if not old_n_fittest == n_fittest:
            diff_out('n_fittest', old_n_fittest, n_fittest)
            io_stat.set_input_common(stat, sec, 'n_fittest', n_fittest)
            logic_change = True
        if not old_slct_func == slct_func:
            diff_out('slct_func', old_slct_func, slct_func)
            io_stat.set_input_common(stat, sec, 'slct_func', slct_func)
            logic_change = True
        if old_slct_func == 'TNM' and slct_func == 'TNM':
            if not old_t_size == t_size:
                diff_out('t_size', old_t_size, t_size)
                io_stat.set_input_common(stat, sec, 't_size', t_size)
                logic_change = True
        elif old_slct_func == 'RLT' and slct_func == 'RLT':
            if not old_a_rlt == a_rlt:
                diff_out('a_rlt', old_a_rlt, a_rlt)
                io_stat.set_input_common(stat, sec, 'a_rlt', a_rlt)
                logic_change = True
            if not old_b_rlt == b_rlt:
                diff_out('b_rlt', old_b_rlt, b_rlt)
                io_stat.set_input_common(stat, sec, 'b_rlt', b_rlt)
                logic_change = True
        if not old_crs_lat == crs_lat:
            diff_out('crs_lat', old_crs_lat, crs_lat)
            io_stat.set_input_common(stat, sec, 'crs_lat', crs_lat)
            logic_change = True
        if not old_nat_diff_tole == nat_diff_tole:
            diff_out('nat_diff_tole', old_nat_diff_tole, nat_diff_tole)
            io_stat.set_input_common(stat, sec, 'nat_diff_tole', nat_diff_tole)
            logic_change = True
        if not old_ntimes == ntimes:
            diff_out('ntimes', old_ntimes, ntimes)
            io_stat.set_input_common(stat, sec, 'ntimes', ntimes)
            logic_change = True
        if not old_sigma_st == sigma_st:
            diff_out('sigma_st', old_sigma_st, sigma_st)
            io_stat.set_input_common(stat, sec, 'sigma_st', sigma_st)
            logic_change = True
        if struc_mode in ['mol', 'mol_bs']:
            if not old_protect_mol_struc == protect_mol_struc:
                diff_out('protect_mol_struc', old_protect_mol_struc, protect_mol_struc)
                io_stat.set_input_common(stat, sec, 'protect_mol_struc', protect_mol_struc)
                logic_change = True
            if not old_rot_max_angle == rot_max_angle:
                diff_out('rot_max_angle', old_rot_max_angle, rot_max_angle)
                io_stat.set_input_common(stat, sec, 'rot_max_angle', rot_max_angle)
                logic_change = True
            if not old_mindist_mol_ea == mindist_mol_ea:
                diff_out('mindist_mol_ea', old_mindist_mol_ea, mindist_mol_ea)
                for i in range(len(mol_file)):
                    io_stat.set_input_common(stat, sec, 'mindist_mol_ea_{}'.format(i+1),
                                             '{}'.format(' '.join(
                                                 str(x) for x in mindist_mol_ea[i])))
        if not old_maxcnt_ea == maxcnt_ea:
            diff_out('maxcnt_ea', old_maxcnt_ea, maxcnt_ea)
            io_stat.set_input_common(stat, sec, 'maxcnt_ea', maxcnt_ea)
            logic_change = True
        if not old_maxgen_ea == maxgen_ea:
            diff_out('maxgen_ea', old_maxgen_ea, maxgen_ea)
            io_stat.set_input_common(stat, sec, 'maxgen_ea', maxgen_ea)
            logic_change = True
        if not old_emax_ea == emax_ea:
            diff_out('emax_ea', old_emax_ea, emax_ea)
            io_stat.set_input_common(stat, sec, 'emax_ea', emax_ea)
            logic_change = True
        if not old_emin_ea == emin_ea:
            diff_out('emin_ea', old_emin_ea, emin_ea)
            io_stat.set_input_common(stat, sec, 'emin_ea', emin_ea)
            logic_change = True
        if algo == 'EA-vc':
            if not old_end_point == end_point:
                raise ValueError('Do not change end_point')

    # ------ VASP
    sec = 'VASP'
    if calc_code == 'VASP':
        if not old_kppvol == kppvol:
            diff_out('kppvol', old_kppvol, kppvol)
            io_stat.set_input_common(stat, sec, 'kppvol', '{}'.format(
                ' '.join(str(x) for x in kppvol)))
            logic_change = True
        if not old_force_gamma == force_gamma:
            diff_out('force_gamma', old_force_gamma, force_gamma)
            io_stat.set_input_common(stat, sec, 'force_gamma', force_gamma)
            logic_change = True

    # ------ QE
    sec = 'QE'
    if calc_code == 'QE':
        if not old_qe_infile == qe_infile:
            raise ValueError('Do not change qe_infile')
        if not old_qe_outfile == qe_outfile:
            raise ValueError('Do not change qe_outfile')
        if not old_kppvol == kppvol:
            diff_out('kppvol', old_kppvol, kppvol)
            io_stat.set_input_common(stat, sec, 'kppvol', '{}'.format(
                ' '.join(str(x) for x in kppvol)))
            logic_change = True
        if not old_force_gamma == force_gamma:
            diff_out('force_gamma', old_force_gamma, force_gamma)
            io_stat.set_input_common(stat, sec, 'force_gamma', force_gamma)
            logic_change = True
        if not old_pv_term == pv_term:
            diff_out('pv_term', old_pv_term, pv_term)
            io_stat.set_input_common(stat, sec, 'pv_term', pv_term)
            logic_change = True

    # ------ OMX
    sec = 'OMX'
    if calc_code == 'OMX':
        if not old_OMX_infile == OMX_infile:
            raise ValueError('Do not change OMX_infile')
        if not old_OMX_outfile == OMX_outfile:
            raise ValueError('Do not change OMX_outfile')
        if not old_kppvol == kppvol:
            diff_out('kppvol', old_kppvol, kppvol)
            io_stat.set_input_common(stat, sec, 'kppvol', '{}'.format(
                ' '.join(str(x) for x in kppvol)))
            logic_change = True
        if not old_force_gamma == force_gamma:
            diff_out('force_gamma', old_force_gamma, force_gamma)
            io_stat.set_input_common(stat, sec, 'force_gamma', force_gamma)
            logic_change = True

    # ------ soiap
    sec = 'soiap'
    if calc_code == 'soiap':
        if not old_soiap_infile == soiap_infile:
            raise ValueError('Do not change soiap_infile')
        if not old_soiap_outfile == soiap_outfile:
            raise ValueError('Do not change soiap_outfile')
        if not old_soiap_cif == soiap_cif:
            raise ValueError('Do not change soiap_cif')

    # ------ lammps
    sec = 'LAMMPS'
    if calc_code == 'LAMMPS':
        if not old_lammps_infile == lammps_infile:
            raise ValueError('Do not change lammps_infile')
        if not old_lammps_outfile == lammps_outfile:
            raise ValueError('Do not change lammps_outfile')
        if not old_lammps_potential == lammps_potential:
            raise ValueError('Do not change lammps_potential')
        if not old_lammps_data == lammps_data:
            raise ValueError('Do not change lammps_data')

    # ------ ASE
    sec = 'ASE'
    if calc_code == 'ASE':
        if not old_ase_python == ase_python:
            raise ValueError('Do not change ase_python')

    # ------ option
    sec = 'option'
    if not old_stop_chkpt == stop_chkpt:
        diff_out('stop_chkpt', old_stop_chkpt, stop_chkpt)
        io_stat.set_input_common(stat, sec, 'stop_chkpt', stop_chkpt)
        logic_change = True
    if not old_load_struc_flag == load_struc_flag:
        diff_out('load_struc_flag', old_load_struc_flag, load_struc_flag)
        io_stat.set_input_common(stat, sec, 'load_struc_flag', load_struc_flag)
        logic_change = True
    if not old_stop_next_struc == stop_next_struc:
        diff_out('stop_next_struc', old_stop_next_struc, stop_next_struc)
        io_stat.set_input_common(stat, sec, 'stop_next_struc', stop_next_struc)
        logic_change = True
    if not old_recalc == recalc:
        diff_out('recalc', old_recalc, recalc)
        io_stat.set_input_common(stat, sec, 'recalc', '{}'.format(
            ' '.join(str(x) for x in recalc)))
        logic_change = True
    if not old_append_struc_ea == append_struc_ea:
        diff_out('append_struc_ea', old_append_struc_ea, append_struc_ea)
        io_stat.set_input_common(stat, sec, 'append_struc_ea', append_struc_ea)
        logic_change = True
    if not old_energy_step_flag == energy_step_flag:
        raise ValueError('Do not change energy_step_flag')
    if not old_struc_step_flag == struc_step_flag:
        raise ValueError('Do not change struc_step_flag')
    if not old_force_step_flag == force_step_flag:
        raise ValueError('Do not change force_step_flag')
    if not old_stress_step_flag == stress_step_flag:
        raise ValueError('Do not change stress_step_flag')

    # ---------- save stat if necessary
    if logic_change:
        io_stat.write_stat(stat)


def diff_out(var_str, old_var, var):
    logger.info(f'Changed {var_str} from {old_var} to {var}')
