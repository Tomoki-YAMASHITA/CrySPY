'''
Read input from cryspy.in
'''

import configparser
import os

from . import io_stat


def readin():
    # ---------- read cryspy.in
    if not os.path.isfile('cryspy.in'):
        raise IOError('Could not find cryspy.in file')
    config = configparser.ConfigParser()
    config.read('cryspy.in')

    # ---------- basic
    # ------ global declaration
    global algo, calc_code, tot_struc, natot
    global atype, nat, nstage, njob, jobcmd, jobfile
    # ------ read intput variables
    algo = config.get('basic', 'algo')
    if algo not in ['RS', 'BO', 'LAQA', 'EA']:
        raise NotImplementedError('algo must be RS, BO, LAQA, or EA')
    calc_code = config.get('basic', 'calc_code')
    if algo == 'LAQA':
        if not calc_code == 'VASP':
            raise NotImplementedError('LAQA: only VASP for now')
    if calc_code not in ['VASP', 'QE', 'soiap', 'LAMMPS']:
        raise NotImplementedError(
            'calc_code must be VASP, QE, soiap, or LAMMPS')
    tot_struc = config.getint('basic', 'tot_struc')
    if tot_struc <= 0:
        raise ValueError('tot_struc <= 0, check tot_struc')
    natot = config.getint('basic', 'natot')
    if natot <= 0:
        raise ValueError('natot <= 0, check natot')
    atype = config.get('basic', 'atype')
    atype = [a for a in atype.split()]    # list
    nat = config.get('basic', 'nat')
    nat = [int(x) for x in nat.split()]    # character --> integer
    if not len(nat) == len(atype):
        raise ValueError('not len(nat) == len(atype), check atype and nat')
    if not sum(nat) == natot:
        raise ValueError('not sum(nat) == natot, check natot and nat')
    nstage = config.getint('basic', 'nstage')
    if nstage <= 0:
        raise ValueError('nstage <= 0, check nstage')
    if algo == 'LAQA':
        if not nstage == 1:
            raise ValueError('nstage shoud be 1 in LAQA')
    njob = config.getint('basic', 'njob')
    if njob <= 0:
        raise ValueError('njob <= 0, check njob')
    jobcmd = config.get('basic', 'jobcmd')
    jobfile = config.get('basic', 'jobfile')

    # ---------- BO
    if algo == 'BO':
        # ------ global declaration
        global nselect_bo, score, num_rand_basis, cdev, dscrpt
        global fp_rmin, fp_rmax, fp_npoints, fp_sigma
        global max_select_bo, manual_select_bo
        # ------ read intput variables
        nselect_bo = config.getint('BO', 'nselect_bo')
        if nselect_bo <= 0:
            raise ValueError('nselect_bo <= 0, check nselect_bo')
        elif tot_struc < nselect_bo:
            raise ValueError('tot_struc < nselect_bo, check nselect_bo')
        score = config.get('BO', 'score')
        if score == 'TS' or score == 'EI' or score == 'PI':
            pass
        else:
            raise ValueError('score must be TS, EI, or PI, check score')
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
            pass
        else:
            raise NotImplementedError('Now FP only')
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
            raise ValueError('fp_rmin < 0, check fp_rmin')
        if fp_rmax < fp_rmin:
            raise ValueError('fp_rmax < fp_rmin, check fp_rmin and fp_rmax')
        try:
            fp_npoints = config.getint('BO', 'fp_npoints')
        except configparser.NoOptionError:
            fp_npoints = 50
        if fp_npoints <= 0:
            raise ValueError('fp_npoints <= 0, check fp_npoints')
        try:
            fp_sigma = config.getfloat('BO', 'fp_sigma')
        except configparser.NoOptionError:
            fp_sigma = 0.2
        if fp_sigma < 0:
            raise ValueError('fp_sigma < 0, check fp_sigma')
        # -- BO option
        try:
            max_select_bo = config.getint('BO', 'max_select_bo')
        except configparser.NoOptionError:
            max_select_bo = 0
        if max_select_bo < 0:
            raise ValueError('max_select_bo must be non-negative int')
        try:
            manual_select_bo = config.get('BO', 'manual_select_bo')
            manual_select_bo = [int(x) for x in manual_select_bo.split()]
        except configparser.NoOptionError:
            manual_select_bo = []
        if manual_select_bo:
            for i in manual_select_bo:
                if not 0 <= i < tot_struc:
                    raise ValueError('manual_select_bo must be'
                                     ' non-negative int'
                                     ' and less than tot_struc')

    # ---------- LAQA
    if algo == 'LAQA':
        # ------ global declaration
        global nselect_laqa, weight_laqa
        # ------ read intput variables
        nselect_laqa = config.getint('LAQA', 'nselect_laqa')
        try:
            weight_laqa = config.getfloat('LAQA', 'weight_laqa')
        except configparser.NoOptionError:
            weight_laqa = 1.0

    # ---------- EA
    # EA part is written below option section

    # ---------- lattice
    # ------ global declaration
    global minlen, maxlen, dangle, mindist
    # ------ read intput variables
    minlen = config.getfloat('lattice', 'minlen')
    maxlen = config.getfloat('lattice', 'maxlen')
    dangle = config.getfloat('lattice', 'dangle')
    if minlen <= 0.0:
        raise ValueError('minlen must be positive')
    if minlen > maxlen:
        raise ValueError('minlen > maxlen')
    if dangle <= 0.0:
        raise ValueError('dangle < 0.0, dangle must be positive')
    mindist = []
    for i in range(len(atype)):
        tmp = config.get('lattice', 'mindist_{}'.format(i+1))
        tmp = [float(x) for x in tmp.split()]    # character --> float
        if not len(tmp) == len(atype):
            raise ValueError('not len(mindist_{}) == len(atype)'.format(i+1))
        mindist.append(tmp)
    # -- check symmetric matrix
    for i in range(len(mindist)):
        for j in range(len(mindist)):
            if i < j:
                if not mindist[i][j] == mindist[j][i]:
                    raise ValueError('mindist is not symmetric. ({}, {}) -->'
                                     ' {}, ({}, {}) --> {}'.format(
                                         i, j, mindist[i][j],
                                         j, i, mindist[j][i]))

    # ---------- global declaration for comman part in calc_code
    global kppvol, kpt_flag, force_gamma

    # ---------- VASP
    if calc_code == 'VASP':
        # ------ read intput variables
        kpt_flag = True
        kppvol = config.get('VASP', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            raise ValueError('not len(kppvol) == nstage,'
                             ' check kppvol and nstage')
        try:
            force_gamma = config.getboolean('VASP', 'force_gamma')
        except configparser.NoOptionError:
            force_gamma = False

    # ---------- QE
    elif calc_code == 'QE':
        # ------ global declaration
        global qe_infile, qe_outfile
        # ------ read intput variables
        kpt_flag = True
        qe_infile = config.get('QE', 'qe_infile')
        qe_outfile = config.get('QE', 'qe_outfile')
        kppvol = config.get('QE', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            raise ValueError('not len(kppvol) == nstage,'
                             ' check kppvol and nstage')
        try:
            force_gamma = config.getboolean('QE', 'force_gamma')
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
    else:
        raise NotImplementedError('calc_code must be VASP, QE, soiap,'
                                  ' or LAMMPS')

    # ---------- option
    # ------ global declaration
    global maxcnt, stop_chkpt, symprec, spgnum
    global load_struc_flag, stop_next_struc, recalc
    global append_struc_ea
    global energy_step_flag, struc_step_flag, fs_step_flag

    # ------ read intput variables
    try:
        maxcnt = config.getint('option', 'maxcnt')
    except (configparser.NoOptionError, configparser.NoSectionError):
        maxcnt = 50
    try:
        stop_chkpt = config.getint('option', 'stop_chkpt')
    except (configparser.NoOptionError, configparser.NoSectionError):
        stop_chkpt = 0
    try:
        symprec = config.getfloat('option', 'symprec')
    except (configparser.NoOptionError, configparser.NoSectionError):
        symprec = 0.001
    try:
        spgnum = config.get('option', 'spgnum')
    except (configparser.NoOptionError, configparser.NoSectionError):
        spgnum = 'all'
    if spgnum == '0':
        spgnum = 0
    elif spgnum == 'all':
        pass
    else:
        spgnum = spglist(spgnum)
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
                raise ValueError('recalc must be non-negative int'
                                 ' and less than tot_struc')
    try:
        append_struc_ea = config.getboolean('option', 'append_struc_ea')
    except (configparser.NoOptionError, configparser.NoSectionError):
        append_struc_ea = False
    try:
        energy_step_flag = config.getboolean('option', 'energy_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            energy_step_flag = False
    except (configparser.NoOptionError, configparser.NoSectionError):
        energy_step_flag = False
    try:
        struc_step_flag = config.getboolean('option', 'struc_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            struc_step_flag = False
    except (configparser.NoOptionError, configparser.NoSectionError):
        struc_step_flag = False
    try:
        fs_step_flag = config.getboolean('option', 'fs_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            fs_step_flag = False
    except (configparser.NoOptionError, configparser.NoSectionError):
        fs_step_flag = False
    if algo == 'LAQA':
        fs_step_flag = True

    # ---------- EA
    if algo == 'EA' or append_struc_ea:
        # ------ global declaration
        global n_pop, n_crsov, n_perm, n_strain, n_rand, n_elite
        global fit_reverse, n_fittest
        global slct_func, t_size, a_rlt, b_rlt
        global crs_lat, crs_func, nat_diff_tole, ntimes, sigma_st,  maxcnt_ea
        global maxgen_ea
        # global restart_gen
        # ------ read intput variables
        # -- number of structures
        n_pop = config.getint('EA', 'n_pop')
        if n_pop <= 0:
            raise ValueError('n_pop must be positive int')
        n_crsov = config.getint('EA', 'n_crsov')
        if n_crsov < 0:
            raise ValueError('n_crsov must be zero or positive int')
        n_perm = config.getint('EA', 'n_perm')
        if n_perm < 0:
            raise ValueError('n_perm must be zero or positive int')
        if n_perm != 0 and len(atype) == 1:
            raise ValueError('When the number of atom type is 1,'
                             ' n_perm must be 0')
        n_strain = config.getint('EA', 'n_strain')
        if n_strain < 0:
            raise ValueError('n_strain must be zero or positive int')
        n_rand = config.getint('EA', 'n_rand')
        if n_rand < 0:
            raise ValueError('n_rand must be zero or positive int')
        if n_crsov + n_perm + n_strain + n_rand != n_pop:
            raise ValueError('n_crsov + n_perm + n_strain + n_rand'
                             ' must be n_pop')
        n_elite = config.getint('EA', 'n_elite')
        if n_elite < 0:
            raise ValueError('n_elite must be non-negative int')
        # -- fittest
        try:
            fit_reverse = config.getboolean('EA', 'fit_reverse')
        except configparser.NoOptionError:
            fit_reverse = False
        try:
            n_fittest = config.getint('EA', 'n_fittest')
        except configparser.NoOptionError:
            n_fittest = 0
        if n_fittest < 0:
            raise ValueError('n_fittest must be zero or positive int')
        # -- select function
        slct_func = config.get('EA', 'slct_func')
        if slct_func not in ['TNM', 'RLT']:
            raise ValueError('slct_func must be TNM or RLT')
        if slct_func == 'TNM':
            try:
                t_size = config.getint('EA', 't_size')
            except configparser.NoOptionError:
                t_size = 3
            if t_size < 2:
                raise ValueError('t_size must be greater than or equal to 2')
        elif slct_func == 'RLT':
            try:
                a_rlt = config.getfloat('EA', 'a_rlt')
            except configparser.NoOptionError:
                a_rlt = 2.0
            try:
                b_rlt = config.getfloat('EA', 'b_rlt')
            except configparser.NoOptionError:
                b_rlt = 1.0
        # -- crossover
        try:
            crs_lat = config.get('EA', 'crs_lat')
        except configparser.NoOptionError:
            crs_lat = 'equal'
        if crs_lat not in ['equal', 'random']:
            raise ValueError('crs_lat must be equal or random')
        try:
            crs_func = config.get('EA', 'crs_func')
        except configparser.NoOptionError:
            crs_func = 'OP'
        if crs_func not in ['OP', 'TP']:
            raise ValueError('crs_func must be OP or TP')
        try:
            nat_diff_tole = config.getint('EA', 'nat_diff_tole')
        except configparser.NoOptionError:
            nat_diff_tole = 4
        if nat_diff_tole < 0:
            raise ValueError('nat_diff_tole must be nen-negative int')
        # -- permutation
        try:
            ntimes = config.getint('EA', 'ntimes')
        except configparser.NoOptionError:
            ntimes = 1
        if ntimes <= 0:
            raise ValueError('ntimes must be positive int')
        try:
            sigma_st = config.getfloat('EA', 'sigma_st')
        except configparser.NoOptionError:
            sigma_st = 0.5
        if sigma_st <= 0:
            raise ValueError('simga_st must be positive float')
        # -- common
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
            raise ValueError('maxgen_ea must be non-negative int')
        # # -- restart option
        # try:
        #     restart_gen = config.getint('EA', 'restart_gen')
        # except configparser.NoOptionError:
        #     restart_gen = 0


def spglist(spgnum):
    tmpspg = []
    for c in spgnum.split():
        if '-' in c:
            if not len(c.split('-')) == 2:
                raise ValueError('Wrong input in spgnum. ')
            istart = int(c.split('-')[0])
            iend = int(c.split('-')[1])+1
            if istart < 0 or 230 < istart:
                raise ValueError('spgnum must be 1 -- 230')
            if iend < 0 or 231 < iend:
                raise ValueError('spgnum must be 1 -- 230')
            for i in range(istart, iend):
                if i not in tmpspg:
                    tmpspg.append(i)
        else:
            if int(c) < 0 or 230 < int(c):
                raise ValueError('spgnum must be 1 -- 230')
            if not int(c) in tmpspg:
                tmpspg += [int(c)]
    return tmpspg


def writeout():
    # ---------- write input data in output file
    print('Write input data in cryspy.out')
    with open('cryspy.out', 'a') as fout:
        fout.write('# ---------- Read cryspy.in (at 1st run)\n')
        fout.write('# ------ basic section\n')
        fout.write('algo = {}\n'.format(algo))
        fout.write('calc_code = {}\n'.format(calc_code))
        fout.write('tot_struc = {}\n'.format(tot_struc))
        fout.write('natot = {}\n'.format(natot))
        fout.write('atype = {}\n'.format(' '.join(a for a in atype)))
        fout.write('nat = {}\n'.format(' '.join(str(b) for b in nat)))
        fout.write('nstage = {}\n'.format(nstage))
        fout.write('njob = {}\n'.format(njob))
        fout.write('jobcmd = {}\n'.format(jobcmd))
        fout.write('jobfile = {}\n'.format(jobfile))

        # ------ BO
        if algo == 'BO':
            fout.write('# ------ BO section\n')
            fout.write('nselect_bo = {}\n'.format(nselect_bo))
            fout.write('score = {}\n'.format(score))
            fout.write('num_rand_basis = {}\n'.format(num_rand_basis))
            fout.write('cdev = {}\n'.format(cdev))
            fout.write('dscrpt = {}\n'.format(dscrpt))
            fout.write('fp_rmin = {}\n'.format(fp_rmin))
            fout.write('fp_rmax = {}\n'.format(fp_rmax))
            fout.write('fp_npoints = {}\n'.format(fp_npoints))
            fout.write('fp_sigma = {}\n'.format(fp_sigma))
            fout.write('max_select_bo = {}\n'.format(max_select_bo))
            fout.write('manual_select_bo = {}\n'.format(
                ' '.join(str(x) for x in manual_select_bo)))

        # ------ LAQA
        if algo == 'LAQA':
            fout.write('# ------ LAQA section\n')
            fout.write('nselect_laqa = {}\n'.format(nselect_laqa))
            fout.write('weight_laqa = {}\n'.format(weight_laqa))

        # ------ EA
        if algo == 'EA' or append_struc_ea:
            fout.write('# ------ EA section\n')
            fout.write('n_pop = {}\n'.format(n_pop))
            fout.write('n_crsov = {}\n'.format(n_crsov))
            fout.write('n_perm = {}\n'.format(n_perm))
            fout.write('n_strain = {}\n'.format(n_strain))
            fout.write('n_rand = {}\n'.format(n_rand))
            fout.write('n_elite = {}\n'.format(n_elite))
            fout.write('fit_reverse = {}\n'.format(fit_reverse))
            fout.write('n_fittest = {}\n'.format(n_fittest))
            fout.write('slct_func = {}\n'.format(slct_func))
            if slct_func == 'TNM':
                fout.write('t_size = {}\n'.format(t_size))
            elif slct_func == 'RLT':
                fout.write('a_rlt = {}\n'.format(a_rlt))
                fout.write('b_rlt = {}\n'.format(b_rlt))
            fout.write('crs_lat = {}\n'.format(crs_lat))
            fout.write('crs_func = {}\n'.format(crs_func))
            fout.write('nat_diff_tole = {}\n'.format(nat_diff_tole))
            fout.write('ntimes = {}\n'.format(ntimes))
            fout.write('sigma_st = {}\n'.format(sigma_st))
            fout.write('maxcnt_ea = {}\n'.format(maxcnt_ea))
            fout.write('maxgen_ea = {}\n'.format(maxgen_ea))
#            fout.write('restart_gen = {}\n'.format(restart_gen))

        # ------ lattice
        fout.write('# ------ lattice section\n')
        fout.write('minlen = {}\n'.format(minlen))
        fout.write('maxlen = {}\n'.format(maxlen))
        fout.write('dangle = {}\n'.format(dangle))
        for i in range(len(atype)):
            fout.write('mindist_{0} = {1}\n'.format(
                i+1, ' '.join(str(c) for c in mindist[i])))

        # ------ VASP
        if calc_code == 'VASP':
            fout.write('# ------ VASP section\n')
            fout.write('kppvol = {}\n'.format(
                ' '.join(str(c) for c in kppvol)))
            fout.write('force_gamma = {}\n'.format(force_gamma))

        # ------- QE
        if calc_code == 'QE':
            fout.write('# ------ QE section\n')
            fout.write('qe_infile = {}\n'.format(qe_infile))
            fout.write('qe_outfile = {}\n'.format(qe_outfile))
            fout.write('kppvol = {}\n'.format(
                ' '.join(str(c) for c in kppvol)))
            fout.write('force_gamma = {}\n'.format(force_gamma))

        # ------ soiap
        if calc_code == 'soiap':
            fout.write('# ------ soiap section\n')
            fout.write('soiap_infile = {}\n'.format(soiap_infile))
            fout.write('soiap_outfile = {}\n'.format(soiap_outfile))
            fout.write('soiap_cif = {}\n'.format(soiap_cif))

        # ------ lammps
        if calc_code == 'LAMMPS':
            fout.write('# ------ lammps section\n')
            fout.write('lammps_infile = {}\n'.format(lammps_infile))
            fout.write('lammps_outfile = {}\n'.format(lammps_outfile))
            fout.write('lammps_potential = {}\n'.format(
                ' '.join(lammps_potential)))
            fout.write('lammps_data = {}\n'.format(lammps_data))

        # ------ option
        fout.write('# ------ option section\n')
        fout.write('maxcnt = {}\n'.format(maxcnt))
        fout.write('stop_chkpt = {}\n'.format(stop_chkpt))
        fout.write('symprec = {}\n'.format(symprec))
        if spgnum == 0 or spgnum == 'all':
            fout.write('spgnum = {}\n'.format(spgnum))
        else:
            fout.write('spgnum = {}\n'.format(
                ' '.join(str(d) for d in spgnum)))
        fout.write('load_struc_flag = {}\n'.format(load_struc_flag))
        fout.write('stop_next_struc = {}\n'.format(stop_next_struc))
        fout.write('recalc = {}\n'.format(' '.join(str(x) for x in recalc)))
        fout.write('append_struc_ea = {}\n'.format(append_struc_ea))
        fout.write('energy_step_flag = {}\n'.format(energy_step_flag))
        fout.write('struc_step_flag = {}\n'.format(struc_step_flag))
        fout.write('fs_step_flag = {}\n'.format(fs_step_flag))
        fout.write('\n\n')


def save_stat(stat):
    print('Save input data in cryspy.stat')
    # ---------- basic
    stat.set('input', 'algo', '{}'.format(algo))
    stat.set('input', 'calc_code', '{}'.format(calc_code))
    stat.set('input', 'tot_struc', '{}'.format(tot_struc))
    stat.set('input', 'natot', '{}'.format(natot))
    stat.set('input', 'atype', '{}'.format(' '.join(a for a in atype)))
    stat.set('input', 'nat', '{}'.format(' '.join(str(b) for b in nat)))
    stat.set('input', 'nstage', '{}'.format(nstage))
    stat.set('input', 'njob', '{}'.format(njob))
    stat.set('input', 'jobcmd', '{}'.format(jobcmd))
    stat.set('input', 'jobfile', '{}'.format(jobfile))

    # ---------- BO
    if algo == 'BO':
        stat.set('input', 'nselect_bo', '{}'.format(nselect_bo))
        stat.set('input', 'score', '{}'.format(score))
        stat.set('input', 'num_rand_basis', '{}'.format(num_rand_basis))
        stat.set('input', 'cdev', '{}'.format(cdev))
        stat.set('input', 'dscrpt', '{}'.format(dscrpt))
        stat.set('input', 'fp_rmin', '{}'.format(fp_rmin))
        stat.set('input', 'fp_rmax', '{}'.format(fp_rmax))
        stat.set('input', 'fp_npoints', '{}'.format(fp_npoints))
        stat.set('input', 'fp_sigma', '{}'.format(fp_sigma))
        stat.set('input', 'max_select_bo', '{}'.format(max_select_bo))
        stat.set('input', 'manual_select_bo', '{}'.format(
            ' '.join(str(x) for x in manual_select_bo)))

    # ---------- LAQA
    if algo == 'LAQA':
        stat.set('input', 'nselect_laqa', '{}'.format(nselect_laqa))
        stat.set('input', 'weight_laqa', '{}'.format(weight_laqa))

    # ---------- EA
    elif algo == 'EA' or append_struc_ea:
        stat.set('input', 'n_pop', '{}'.format(n_pop))
        stat.set('input', 'n_crsov', '{}'.format(n_crsov))
        stat.set('input', 'n_perm', '{}'.format(n_perm))
        stat.set('input', 'n_strain', '{}'.format(n_strain))
        stat.set('input', 'n_rand', '{}'.format(n_rand))
        stat.set('input', 'n_elite', '{}'.format(n_elite))
        stat.set('input', 'fit_reverse', '{}'.format(fit_reverse))
        stat.set('input', 'n_fittest', '{}'.format(n_fittest))
        stat.set('input', 'slct_func', '{}'.format(slct_func))
        if slct_func == 'TNM':
            stat.set('input', 't_size', '{}'.format(t_size))
        elif slct_func == 'RLT':
            stat.set('input', 'a_rlt', '{}'.format(a_rlt))
            stat.set('input', 'b_rlt', '{}'.format(b_rlt))
        stat.set('input', 'crs_func', '{}'.format(crs_func))
        stat.set('input', 'crs_lat', '{}'.format(crs_lat))
        stat.set('input', 'nat_diff_tole', '{}'.format(nat_diff_tole))
        stat.set('input', 'ntimes', '{}'.format(ntimes))
        stat.set('input', 'sigma_st', '{}'.format(sigma_st))
        stat.set('input', 'maxcnt_ea', '{}'.format(maxcnt_ea))
        stat.set('input', 'maxgen_ea', '{}'.format(maxgen_ea))

    # ---------- lattice
    stat.set('input', 'minlen', '{}'.format(minlen))
    stat.set('input', 'maxlen', '{}'.format(maxlen))
    stat.set('input', 'dangle', '{}'.format(dangle))
    for i in range(len(atype)):
        stat.set('input', 'mindist_{}'.format(i+1),
                 '{}'.format(' '.join(str(c) for c in mindist[i])))

    # ---------- VASP
    if calc_code == 'VASP':
        stat.set('input', 'kppvol',
                 '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('input', 'force_gamma', '{}'.format(force_gamma))

    # ---------- QE
    if calc_code == 'QE':
        stat.set('input', 'qe_infile', '{}'.format(qe_infile))
        stat.set('input', 'qe_outfile', '{}'.format(qe_outfile))
        stat.set('input', 'kppvol',
                 '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('input', 'force_gamma', '{}'.format(force_gamma))

    # ---------- soiap
    if calc_code == 'soiap':
        stat.set('input', 'soiap_infile', '{}'.format(soiap_infile))
        stat.set('input', 'soiap_outfile', '{}'.format(soiap_outfile))
        stat.set('input', 'soiap_cif', '{}'.format(soiap_cif))

    # ---------- lammps
    if calc_code == 'LAMMPS':
        stat.set('input', 'lammps_infile', '{}'.format(lammps_infile))
        stat.set('input', 'lammps_outfile', '{}'.format(lammps_outfile))
        stat.set('input', 'lammps_potential',
                 '{}'.format(' '.join(lammps_potential)))
        stat.set('input', 'lammps_data', '{}'.format(lammps_data))

    # ---------- option
    stat.set('input', 'maxcnt', '{}'.format(maxcnt))
    stat.set('input', 'stop_chkpt', '{}'.format(stop_chkpt))
    stat.set('input', 'symprec', '{}'.format(symprec))
    if spgnum == 0 or spgnum == 'all':
        stat.set('input', 'spgnum', '{}'.format(spgnum))
    else:
        stat.set('input', 'spgnum',
                 '{}'.format(' '.join(str(d) for d in spgnum)))
    stat.set('input', 'load_struc_flag', '{}'.format(load_struc_flag))
    stat.set('input', 'stop_next_struc', '{}'.format(stop_next_struc))
    stat.set('input', 'recalc', '{}'.format(' '.join(str(x) for x in recalc)))
    stat.set('input', 'append_struc_ea', '{}'.format(append_struc_ea))
    stat.set('input', 'energy_step_flag', '{}'.format(energy_step_flag))
    stat.set('input', 'struc_step_flag', '{}'.format(struc_step_flag))
    stat.set('input', 'fs_step_flag', '{}'.format(fs_step_flag))

    # ---------- write stat
    io_stat.write_stat(stat)


def diffinstat(stat):
    logic_change = False

    # ---------- old input
    # ------ basic
    old_algo = stat.get('input', 'algo')
    old_calc_code = stat.get('input', 'calc_code')
    old_tot_struc = stat.getint('input', 'tot_struc')
    old_natot = stat.getint('input', 'natot')
    old_atype = stat.get('input', 'atype')
    old_atype = [a for a in old_atype.split()]    # list
    old_nat = stat.get('input', 'nat')
    old_nat = [int(x) for x in old_nat.split()]    # str --> int list
    old_nstage = stat.getint('input', 'nstage')
    old_njob = stat.getint('input', 'njob')
    old_jobcmd = stat.get('input', 'jobcmd')
    old_jobfile = stat.get('input', 'jobfile')

    # ------ BO
    if old_algo == 'BO':
        old_nselect_bo = stat.getint('input', 'nselect_bo')
        old_score = stat.get('input', 'score')
        old_num_rand_basis = stat.getint('input', 'num_rand_basis')
        old_cdev = stat.getfloat('input', 'cdev')
        old_dscrpt = stat.get('input', 'dscrpt')
        old_fp_rmin = stat.getfloat('input', 'fp_rmin')
        old_fp_rmax = stat.getfloat('input', 'fp_rmax')
        old_fp_npoints = stat.getint('input', 'fp_npoints')
        old_fp_sigma = stat.getfloat('input', 'fp_sigma')
        old_max_select_bo = stat.getint('input', 'max_select_bo')
        old_manual_select_bo = stat.get('input', 'manual_select_bo')
        old_manual_select_bo = [int(x) for x in old_manual_select_bo.split()]

    # ------ LAQA
    if old_algo == 'LAQA':
        old_nselect_laqa = stat.getint('input', 'nselect_laqa')
        old_weight_laqa = stat.getfloat('input', 'weight_laqa')

    # ------ EA
    if old_algo == 'EA':
        old_n_pop = stat.getint('input', 'n_pop')
        old_n_crsov = stat.getint('input', 'n_crsov')
        old_n_perm = stat.getint('input', 'n_perm')
        old_n_strain = stat.getint('input', 'n_strain')
        old_n_rand = stat.getint('input', 'n_rand')
        old_n_elite = stat.getint('input', 'n_elite')
        old_fit_reverse = stat.getboolean('input', 'fit_reverse')
        old_n_fittest = stat.getint('input', 'n_fittest')
        old_slct_func = stat.get('input', 'slct_func')
        if old_slct_func == 'TNM':
            old_t_size = stat.getint('input', 't_size')
        elif old_slct_func == 'RLT':
            old_a_rlt = stat.getfloat('input', 'a_rlt')
            old_b_rlt = stat.getfloat('input', 'b_rlt')
        old_crs_lat = stat.get('input', 'crs_lat')
        old_crs_func = stat.get('input', 'crs_func')
        old_nat_diff_tole = stat.getint('input', 'nat_diff_tole')
        old_ntimes = stat.getint('input', 'ntimes')
        old_sigma_st = stat.getfloat('input', 'sigma_st')
        old_maxcnt_ea = stat.getint('input', 'maxcnt_ea')
        old_maxgen_ea = stat.getint('input', 'maxgen_ea')
        # old_restart_gen = stat.get('input', 'restart_gen')

    # ------ lattice
    old_minlen = stat.getfloat('input', 'minlen')
    old_maxlen = stat.getfloat('input', 'maxlen')
    old_dangle = stat.getfloat('input', 'dangle')
    old_mindist = []
    for i in range(len(atype)):
        tmp = stat.get('input', 'mindist_{}'.format(i+1))
        tmp = [float(x) for x in tmp.split()]    # character --> float
        old_mindist.append(tmp)

    # ------ VASP
    if old_calc_code == 'VASP':
        old_kppvol = stat.get('input', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]    # int list
        old_force_gamma = stat.getboolean('input', 'force_gamma')

    # ------ QE
    if old_calc_code == 'QE':
        old_qe_infile = stat.get('input', 'qe_infile')
        old_qe_outfile = stat.get('input', 'qe_outfile')
        old_kppvol = stat.get('input', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]  # int list
        old_force_gamma = stat.getboolean('input', 'force_gamma')

    # ------ soiap
    if old_calc_code == 'soiap':
        old_soiap_infile = stat.get('input', 'soiap_infile')
        old_soiap_outfile = stat.get('input', 'soiap_outfile')
        old_soiap_cif = stat.get('input', 'soiap_cif')

    # ------ lammps
    if old_calc_code == 'LAMMPS':
        old_lammps_infile = stat.get('input', 'lammps_infile')
        old_lammps_outfile = stat.get('input', 'lammps_outfile')
        old_lammps_potential = stat.get('input', 'lammps_potential')
        old_lammps_potential = old_lammps_potential.split()    # str --> list
        if old_lammps_potential == 'None':    # 'None' is just character here
            old_lammps_potential = None
        old_lammps_data = stat.get('input', 'lammps_data')

    # ------ option
    old_maxcnt = stat.getint('input', 'maxcnt')
    old_stop_chkpt = stat.getint('input', 'stop_chkpt')
    old_symprec = stat.getfloat('input', 'symprec')
    old_spgnum = stat.get('input', 'spgnum')
    if old_spgnum == '0':
        old_spgnum = 0
    elif not old_spgnum == 'all':
        old_spgnum = [int(x) for x in old_spgnum.split()]    # int list
    old_load_struc_flag = stat.getboolean('input', 'load_struc_flag')
    old_stop_next_struc = stat.getboolean('input', 'stop_next_struc')
    old_recalc = stat.get('input', 'recalc')
    old_recalc = [int(x) for x in old_recalc.split()]    # int list
    old_append_struc_ea = stat.getboolean('input', 'append_struc_ea')
    old_energy_step_flag = stat.getboolean('input', 'energy_step_flag')
    old_struc_step_flag = stat.getboolean('input', 'struc_step_flag')
    old_fs_step_flag = stat.getboolean('input', 'fs_step_flag')

    # ---------- check difference
    # ------ basic
    if not old_algo == algo:
        raise ValueError('Do not change algo')
    if not old_calc_code == calc_code:
        raise ValueError('Do not change calc code')
    if not old_tot_struc == tot_struc:
        if algo == 'EA':
            raise ValueError('Do not change tot_struc in EA')
        diff_out('tot_struc', old_tot_struc, tot_struc)
        io_stat.set_input_common(stat, 'tot_struc', tot_struc)
        logic_change = True
    if not old_natot == natot:
        raise ValueError('Do not change natot')
    if not old_atype == atype:
        raise ValueError('Do not change atype')
    if not old_nat == nat:
        raise ValueError('Do not change nat')
    if not old_nstage == nstage:
        diff_out('nstage', old_nstage, nstage)
        io_stat.set_input_common(stat, 'nstage', nstage)
        logic_change = True
    if not old_njob == njob:
        diff_out('njob', old_njob, njob)
        io_stat.set_input_common(stat, 'njob', njob)
        logic_change = True
    if not old_jobcmd == jobcmd:
        diff_out('jobcmd', old_jobcmd, jobcmd)
        io_stat.set_input_common(stat, 'jobcmd', jobcmd)
        logic_change = True
    if not old_jobfile == jobfile:
        diff_out('jobfile', old_jobfile, jobfile)
        io_stat.set_input_common(stat, 'jobfile', jobfile)
        logic_change = True

    # ------ BO
    if algo == 'BO':
        if not old_nselect_bo == nselect_bo:
            diff_out('nselect_bo', old_nselect_bo, nselect_bo)
            io_stat.set_input_common(stat, 'nselect_bo', nselect_bo)
            logic_change = True
        if not old_score == score:
            diff_out('score', old_score, score)
            io_stat.set_input_common(stat, 'score', score)
            logic_change = True
        if not old_num_rand_basis == num_rand_basis:
            diff_out('num_rand_basis', old_num_rand_basis, num_rand_basis)
            io_stat.set_input_common(stat, 'num_rand_basis', num_rand_basis)
            logic_change = True
        if not old_cdev == cdev:
            diff_out('cdev', old_cdev, cdev)
            io_stat.set_input_common(stat, 'cdev', cdev)
            logic_change = True
        if not old_dscrpt == dscrpt:
            raise ValueError('Do not change dscrpt')
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
            io_stat.set_input_common(stat, 'max_select_bo', max_select_bo)
            logic_change = True
        if not old_manual_select_bo == manual_select_bo:
            diff_out('manual_select_bo', old_manual_select_bo,
                     manual_select_bo)
            io_stat.set_input_common(stat, 'manual_select_bo',
                                     '{}'.format(' '.join(
                                         str(x) for x in manual_select_bo)))
            logic_change = True

    # ------ LAQA
    if algo == 'LAQA':
        if not old_nselect_laqa == nselect_laqa:
            diff_out('nselect_laqa', old_nselect_laqa, nselect_laqa)
            io_stat.set_input_common(stat, 'nselect_laqa', nselect_laqa)
            logic_change = True
        if not old_weight_laqa == weight_laqa:
            diff_out('weight_laqa', old_weight_laqa, weight_laqa)
            io_stat.set_input_common(stat, 'weight_laqa', weight_laqa)
            logic_change = True

    # ------ EA
    if algo == 'EA':
        if not old_n_pop == n_pop:
            diff_out('n_pop', old_n_pop, n_pop)
            io_stat.set_input_common(stat, 'n_pop', n_pop)
            logic_change = True
        if not old_n_crsov == n_crsov:
            diff_out('n_crsov', old_n_crsov, n_crsov)
            io_stat.set_input_common(stat, 'n_crsov', n_crsov)
            logic_change = True
        if not old_n_perm == n_perm:
            diff_out('n_perm', old_n_perm, n_perm)
            io_stat.set_input_common(stat, 'n_perm', n_perm)
            logic_change = True
        if not old_n_strain == n_strain:
            diff_out('n_strain', old_n_strain, n_strain)
            io_stat.set_input_common(stat, 'n_strain', n_strain)
            logic_change = True
        if not old_n_rand == n_rand:
            diff_out('n_rand', old_n_rand, n_rand)
            io_stat.set_input_common(stat, 'n_rand', n_rand)
            logic_change = True
        if not old_n_elite == n_elite:
            diff_out('n_elite', old_n_elite, n_elite)
            io_stat.set_input_common(stat, 'n_elite', n_elite)
            logic_change = True
        if not old_fit_reverse == fit_reverse:
            raise ValueError('Do not change fit_reverse')
        if not old_n_fittest == n_fittest:
            diff_out('n_fittest', old_n_fittest, n_fittest)
            io_stat.set_input_common(stat, 'n_fittest', n_fittest)
            logic_change = True
        if not old_slct_func == slct_func:
            diff_out('slct_func', old_slct_func, slct_func)
            io_stat.set_input_common(stat, 'slct_func', slct_func)
            logic_change = True
        if old_slct_func == 'TNM' and slct_func == 'TNM':
            if not old_t_size == t_size:
                diff_out('t_size', old_t_size, t_size)
                io_stat.set_input_common(stat, 't_size', t_size)
                logic_change = True
        elif old_slct_func == 'RLT' and slct_func == 'RLT':
            if not old_a_rlt == a_rlt:
                diff_out('a_rlt', old_a_rlt, a_rlt)
                io_stat.set_input_common(stat, 'a_rlt', a_rlt)
                logic_change = True
            if not old_b_rlt == b_rlt:
                diff_out('b_rlt', old_b_rlt, b_rlt)
                io_stat.set_input_common(stat, 'b_rlt', b_rlt)
                logic_change = True
        if not old_crs_func == crs_func:
            diff_out('crs_func', old_crs_func, crs_func)
            io_stat.set_input_common(stat, 'crs_func', crs_func)
            logic_change = True
        if not old_crs_lat == crs_lat:
            diff_out('crs_lat', old_crs_lat, crs_lat)
            io_stat.set_input_common(stat, 'crs_lat', crs_lat)
            logic_change = True
        if not old_nat_diff_tole == nat_diff_tole:
            diff_out('nat_diff_tole', old_nat_diff_tole, nat_diff_tole)
            io_stat.set_input_common(stat, 'nat_diff_tole', nat_diff_tole)
            logic_change = True
        if not old_ntimes == ntimes:
            diff_out('ntimes', old_ntimes, ntimes)
            io_stat.set_input_common(stat, 'ntimes', ntimes)
            logic_change = True
        if not old_sigma_st == sigma_st:
            diff_out('sigma_st', old_sigma_st, sigma_st)
            io_stat.set_input_common(stat, 'sigma_st', sigma_st)
            logic_change = True
        if not old_maxcnt_ea == maxcnt_ea:
            diff_out('maxcnt_ea', old_maxcnt_ea, maxcnt_ea)
            io_stat.set_input_common(stat, 'maxcnt_ea', maxcnt_ea)
            logic_change = True
        if not old_maxgen_ea == maxgen_ea:
            diff_out('maxgen_ea', old_maxgen_ea, maxgen_ea)
            io_stat.set_input_common(stat, 'maxgen_ea', maxgen_ea)
            logic_change = True

    # ------ lattice
    if not old_minlen == minlen:
        diff_out('minlen', old_minlen, minlen)
        io_stat.set_input_common(stat, 'minlen', minlen)
        logic_change = True
    if not old_maxlen == maxlen:
        diff_out('maxlen', old_maxlen, maxlen)
        io_stat.set_input_common(stat, 'maxlen', maxlen)
        logic_change = True
    if not old_dangle == dangle:
        diff_out('dangle', old_dangle, dangle)
        io_stat.set_input_common(stat, 'dangle', dangle)
        logic_change = True
    if not old_mindist == mindist:
        diff_out('mindist', old_mindist, mindist)
        io_stat.set_input_common(stat, 'mindist', mindist)

        for i in range(len(atype)):
            io_stat.set_input_common(stat, 'mindist_{}'.format(i+1),
                                     '{}'.format(' '.join(
                                         str(x) for x in mindist[i])))
        logic_change = True

    # ------ VASP
    if calc_code == 'VASP':
        if not old_kppvol == kppvol:
            diff_out('kppvol', old_kppvol, kppvol)
            io_stat.set_input_common(stat, 'kppvol', '{}'.format(
                ' '.join(str(x) for x in kppvol)))
            logic_change = True
        if not old_force_gamma == force_gamma:
            diff_out('force_gamma', old_force_gamma, force_gamma)
            io_stat.set_input_common(stat, 'force_gamma', force_gamma)
            logic_change = True

    # ------ QE
    if calc_code == 'QE':
        if not old_qe_infile == qe_infile:
            raise ValueError('Do not change qe_infile')
        if not old_qe_outfile == qe_outfile:
            raise ValueError('Do not change qe_outfile')
        if not old_kppvol == kppvol:
            diff_out('kppvol', old_kppvol, kppvol)
            io_stat.set_input_common(stat, 'kppvol', '{}'.format(
                ' '.join(str(x) for x in kppvol)))
            logic_change = True
        if not old_force_gamma == force_gamma:
            diff_out('force_gamma', old_force_gamma, force_gamma)
            io_stat.set_input_common(stat, 'force_gamma', force_gamma)
            logic_change = True

    # ------ soiap
    if calc_code == 'soiap':
        if not old_soiap_infile == soiap_infile:
            raise ValueError('Do not change soiap_infile')
        if not old_soiap_outfile == soiap_outfile:
            raise ValueError('Do not change soiap_outfile')
        if not old_soiap_cif == soiap_cif:
            raise ValueError('Do not change soiap_cif')

    # ------ lammps
    if calc_code == 'LAMMPS':
        if not old_lammps_infile == lammps_infile:
            raise ValueError('Do not change lammps_infile')
        if not old_lammps_outfile == lammps_outfile:
            raise ValueError('Do not change lammps_outfile')
        if not old_lammps_potential == lammps_potential:
            raise ValueError('Do not change lammps_potential')
        if not old_lammps_data == lammps_data:
            raise ValueError('Do not change lammps_data')

    # ------ option
    if not old_maxcnt == maxcnt:
        diff_out('maxcnt', old_maxcnt, maxcnt)
        io_stat.set_input_common(stat, 'maxcnt', maxcnt)
        logic_change = True
    if not old_stop_chkpt == stop_chkpt:
        diff_out('stop_chkpt', old_stop_chkpt, stop_chkpt)
        io_stat.set_input_common(stat, 'stop_chkpt', stop_chkpt)
        logic_change = True
    if not old_symprec == symprec:
        diff_out('symprec', old_symprec, symprec)
        io_stat.set_input_common(stat, 'symprec', symprec)
        logic_change = True
    if not old_spgnum == spgnum:
        diff_out('spgnum', old_spgnum, spgnum)
        if spgnum == 0 or spgnum == 'all':
            io_stat.set_input_common(stat, 'spgnum', spgnum)
        else:
            io_stat.set_input_common(stat, 'spgnum', '{}'.format(
                ' '.join(str(x) for x in spgnum)))
        logic_change = True
    if not old_load_struc_flag == load_struc_flag:
        diff_out('load_struc_flag', old_load_struc_flag, load_struc_flag)
        io_stat.set_input_common(stat, 'load_struc_flag', load_struc_flag)
        logic_change = True
    if not old_stop_next_struc == stop_next_struc:
        diff_out('stop_next_struc', old_stop_next_struc, stop_next_struc)
        io_stat.set_input_common(stat, 'stop_next_struc', stop_next_struc)
        logic_change = True
    if not old_recalc == recalc:
        diff_out('recalc', old_recalc, recalc)
        io_stat.set_input_common(stat, 'recalc', '{}'.format(
            ' '.join(str(x) for x in recalc)))
        logic_change = True
    if not old_append_struc_ea == append_struc_ea:
        diff_out('append_struc_ea', old_append_struc_ea, append_struc_ea)
        io_stat.set_input_common(stat, 'append_struc_ea', append_struc_ea)
        logic_change = True
    if not old_energy_step_flag == energy_step_flag:
        raise ValueError('Do not change energy_step_flag')
    if not old_struc_step_flag == struc_step_flag:
        raise ValueError('Do not change struc_step_flag')
    if not old_fs_step_flag == fs_step_flag:
        raise ValueError('Do not change fs_step_flag')

    # ---------- save stat if necessary
    if logic_change:
        io_stat.write_stat(stat)


def diff_out(var_str, old_var, var):
    print('Changed {0} from {1} to {2}'.format(var_str, old_var, var))
    with open('cryspy.out', 'a') as fout:
        fout.write('\n#### Changed {0} from {1} to {2}\n'.format(
            var_str, old_var, var))
