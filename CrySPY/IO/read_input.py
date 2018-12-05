#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division

import ConfigParser
import os


def readin():
    # ---------- read cryspy.in
    if not os.path.isfile('cryspy.in'):
        raise IOError('Could not find cryspy.in file')
    config = ConfigParser.ConfigParser()
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
        raise NotImplementedError('calc_code must be VASP, QE, soiap, or LAMMPS')
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
        global interval, score, num_rand_basis, cdev, dscrpt
        global fp_rmin, fp_rmax, fp_npoints, fp_sigma
        global maxgen
        # ------ read intput variables
        interval = config.getint('BO', 'interval')
        if interval <= 0:
            raise ValueError('interval <= 0, check interval')
        elif tot_struc < interval:
            raise ValueError('tot_struc < interval, check interval')
        score = config.get('BO', 'score')
        if score == 'TS' or score == 'EI' or score == 'PI':
            pass
        else:
            raise ValueError('score must be TS, EI, or PI, check score')
        try:
            num_rand_basis = config.getint('BO', 'num_rand_basis')
        except ConfigParser.NoOptionError:
            num_rand_basis = 0
        try:
            cdev = config.getfloat('BO', 'cdev')
        except ConfigParser.NoOptionError:
            cdev = 0.001
        dscrpt = config.get('BO', 'dscrpt')
        if dscrpt == 'FP':
            pass
        else:
            raise NotImplementedError('Now FP only')
        # -- parameters for f-fingerprint (optional)
        try:
            fp_rmin = config.getfloat('BO', 'fp_rmin')
        except ConfigParser.NoOptionError:
            fp_rmin = 0.5
        try:
            fp_rmax = config.getfloat('BO', 'fp_rmax')
        except ConfigParser.NoOptionError:
            fp_rmax = 5.0
        if fp_rmin < 0.0:
            raise ValueError('fp_rmin < 0, check fp_rmin')
        if fp_rmax < fp_rmin:
            raise ValueError('fp_rmax < fp_rmin, check fp_rmin and fp_rmax')
        try:
            fp_npoints = config.getint('BO', 'fp_npoints')
        except ConfigParser.NoOptionError:
            fp_npoints = 50
        if fp_npoints <= 0:
            raise ValueError('fp_npoints <= 0, check fp_npoints')
        try:
            fp_sigma = config.getfloat('BO', 'fp_sigma')
        except ConfigParser.NoOptionError:
            fp_sigma = 0.2
        if fp_sigma < 0:
            raise ValueError('fp_sigma < 0, check fp_sigma')
        # -- BO option
        try:
            maxgen = config.getint('BO', 'maxgen')
        except ConfigParser.NoOptionError:
            maxgen = 0
        if maxgen < 0:
            raise ValueError('maxgen must be non-negative int')

    # ---------- LAQA
    if algo == 'LAQA':
        # ------ global declaration
        global nselect, weight_laqa
        # ------ read intput variables
        nselect = config.getint('LAQA', 'nselect')
        try:
            weight_laqa = config.getfloat('LAQA', 'weight_laqa')
        except ConfigParser.NoOptionError:
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
    if dangle < 0.0:
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
                    raise ValueError('mindist is not symmetric. ({}, {}) --> {}, ({}, {}) --> {}'.format(
                        i, j, mindist[i][j], j, i, mindist[j][i]))

    # ---------- global declaration for comman part in calc_code
    global kppvol, kpt_flag, force_gamma

    # ---------- VASP
    if calc_code == 'VASP':
        # ------ read intput variables
        kpt_flag = True
        kppvol = config.get('VASP', 'kppvol')
        kppvol = [int(x) for x in kppvol.split()]    # character --> int
        if not len(kppvol) == nstage:
            raise ValueError('not len(kppvol) == nstage, check kppvol and nstage')
        try:
            force_gamma = config.getboolean('VASP', 'force_gamma')
        except ConfigParser.NoOptionError:
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
            raise ValueError('not len(kppvol) == nstage, check kppvol and nstage')
        try:
            force_gamma = config.getboolean('QE', 'force_gamma')
        except ConfigParser.NoOptionError:
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
        except ConfigParser.NoOptionError:
            lammps_potential = None
        lammps_data = config.get('LAMMPS', 'lammps_data')
        kpt_flag = False
        force_gamma = False
    else:
        raise NotImplementedError('calc_code must be VASP, QE, soiap, or LAMMPS')

    # ---------- option
    # ------ global declaration
    global maxcnt, stop_chkpt, symprec, spgnum
    global load_struc_flag, stop_next_struc, append_struc_ea
    global energy_step_flag, struc_step_flag, fs_step_flag

    # ------ read intput variables
    try:
        maxcnt = config.getint('option', 'maxcnt')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        maxcnt = 200
    try:
        stop_chkpt = config.getint('option', 'stop_chkpt')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        stop_chkpt = 0
    try:
        symprec = config.getfloat('option', 'symprec')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        symprec = 0.001
    try:
        spgnum = config.get('option', 'spgnum')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        spgnum = 'all'
    if spgnum == '0':
        spgnum = 0
    elif spgnum == 'all':
        pass
    else:
        spgnum = spglist(spgnum)
    try:
        load_struc_flag = config.getboolean('option', 'load_struc_flag')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        load_struc_flag = False
    try:
        stop_next_struc = config.getboolean('option', 'stop_next_struc')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        stop_next_struc = False
    try:
        append_struc_ea = config.getboolean('option', 'append_struc_ea')
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        append_struc_ea = False
    try:
        energy_step_flag = config.getboolean('option', 'energy_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            energy_step_flag = False
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        energy_step_flag = False
    try:
        struc_step_flag = config.getboolean('option', 'struc_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            struc_step_flag = False
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
        struc_step_flag = False
    try:
        fs_step_flag = config.getboolean('option', 'fs_step_flag')
        # -- only VASP or QE for now
        if calc_code in ['soiap', 'LAMMPS']:
            fs_step_flag = False
    except (ConfigParser.NoOptionError, ConfigParser.NoSectionError):
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
        global maxgen
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
            raise ValueError('When the number of atom type is 1, n_perm must be 0')
        n_strain = config.getint('EA', 'n_strain')
        if n_strain < 0:
            raise ValueError('n_strain must be zero or positive int')
        n_rand = config.getint('EA', 'n_rand')
        if n_rand < 0:
            raise ValueError('n_rand must be zero or positive int')
        if n_crsov + n_perm + n_strain + n_rand != n_pop:
            raise ValueError('n_crsov + n_perm + n_strain + n_rand must be n_pop')
        n_elite = config.getint('EA', 'n_elite')
        if n_elite < 0:
            raise ValueError('n_elite must be non-negative int')
        # -- fittest
        try:
            fit_reverse = config.getboolean('EA', 'fit_reverse')
        except ConfigParser.NoOptionError:
            fit_reverse = False
        try:
            n_fittest = config.getint('EA', 'n_fittest')
        except ConfigParser.NoOptionError:
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
            except ConfigParser.NoOptionError:
                t_size = 3
            if t_size < 2:
                raise ValueError('t_size must be greater than or equal to 2')
        elif slct_func == 'RLT':
            try:
                a_rlt = config.getfloat('EA', 'a_rlt')
            except ConfigParser.NoOptionError:
                a_rlt = 2.0
            try:
                b_rlt = config.getfloat('EA', 'b_rlt')
            except ConfigParser.NoOptionError:
                b_rlt = 1.0
        # -- crossover
        try:
            crs_lat = config.get('EA', 'crs_lat')
        except ConfigParser.NoOptionError:
            crs_lat = 'equal'
        if crs_lat not in ['equal', 'random']:
            raise ValueError('crs_lat must be equal or random')
        try:
            crs_func = config.get('EA', 'crs_func')
        except ConfigParser.NoOptionError:
            crs_func = 'OP'
        if crs_func not in ['OP', 'TP']:
            raise ValueError('crs_func must be OP or TP')
        try:
            nat_diff_tole = config.getint('EA', 'nat_diff_tole')
        except ConfigParser.NoOptionError:
            nat_diff_tole = 4
        if nat_diff_tole < 0:
            raise ValueError('nat_diff_tole must be nen-negative int')
        # -- permutation
        try:
            ntimes = config.getint('EA', 'ntimes')
        except ConfigParser.NoOptionError:
            ntimes = 1
        if ntimes <= 0:
            raise ValueError('ntimes must be positive int')
        try:
            sigma_st = config.getfloat('EA', 'sigma_st')
        except ConfigParser.NoOptionError:
            sigma_st = 0.5
        if sigma_st <= 0:
            raise ValueError('simga_st must be positive float')
        # -- common
        try:
            maxcnt_ea = config.getint('EA', 'maxcnt_ea')
        except ConfigParser.NoOptionError:
            maxcnt_ea = 100
        # -- EA option
        try:
            maxgen = config.getint('EA', 'maxgen')
        except ConfigParser.NoOptionError:
            maxgen = 0
        if maxgen < 0:
            raise ValueError('maxgen must be non-negative int')
        # # -- restart option
        # try:
        #     restart_gen = config.getint('EA', 'restart_gen')
        # except ConfigParser.NoOptionError:
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
                if not i in tmpspg:
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
            fout.write('interval = {}\n'.format(interval))
            fout.write('score = {}\n'.format(score))
            fout.write('num_rand_basis = {}\n'.format(num_rand_basis))
            fout.write('cdev = {}\n'.format(cdev))
            fout.write('dscrpt = {}\n'.format(dscrpt))
            fout.write('fp_rmin = {}\n'.format(fp_rmin))
            fout.write('fp_rmax = {}\n'.format(fp_rmax))
            fout.write('fp_npoints = {}\n'.format(fp_npoints))
            fout.write('fp_sigma = {}\n'.format(fp_sigma))
            fout.write('maxgen = {}\n'.format(maxgen))

        # ------ LAQA
        if algo == 'LAQA':
            fout.write('# ------ LAQA section\n')
            fout.write('nselect = {}\n'.format(nselect))
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
            fout.write('maxgen = {}\n'.format(maxgen))
#            fout.write('restart_gen = {}\n'.format(restart_gen))

        # ------ lattice
        fout.write('# ------ lattice section\n')
        fout.write('minlen = {}\n'.format(minlen))
        fout.write('maxlen = {}\n'.format(maxlen))
        fout.write('dangle = {}\n'.format(dangle))
        for i in range(len(atype)):
            fout.write('mindist_{0} = {1}\n'.format(i+1, ' '.join(str(c) for c in mindist[i])))

        # ------ VASP
        if calc_code == 'VASP':
            fout.write('# ------ VASP section\n')
            fout.write('kppvol = {}\n'.format(' '.join(str(c) for c in kppvol)))
            fout.write('force_gamma = {}\n'.format(force_gamma))

        # ------- QE
        if calc_code == 'QE':
            fout.write('# ------ QE section\n')
            fout.write('qe_infile = {}\n'.format(qe_infile))
            fout.write('qe_outfile = {}\n'.format(qe_outfile))
            fout.write('kppvol = {}\n'.format(' '.join(str(c) for c in kppvol)))
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
            fout.write('lammps_potential = {}\n'.format(lammps_potential))
            fout.write('lammps_data = {}\n'.format(lammps_data))

        # ------ option
        fout.write('# ------ option section\n')
        fout.write('maxcnt = {}\n'.format(maxcnt))
        fout.write('stop_chkpt = {}\n'.format(stop_chkpt))
        fout.write('symprec = {}\n'.format(symprec))
        if spgnum == 0 or spgnum == 'all':
            fout.write('spgnum = {}\n'.format(spgnum))
        else:
            fout.write('spgnum = {}\n'.format(' '.join(str(d) for d in spgnum)))
        fout.write('load_struc_flag = {}\n'.format(load_struc_flag))
        fout.write('stop_next_struc = {}\n'.format(stop_next_struc))
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
        stat.set('input', 'interval', '{}'.format(interval))
        stat.set('input', 'score', '{}'.format(score))
        stat.set('input', 'num_rand_basis', '{}'.format(num_rand_basis))
        stat.set('input', 'cdev', '{}'.format(cdev))
        stat.set('input', 'dscrpt', '{}'.format(dscrpt))
        stat.set('input', 'fp_rmin', '{}'.format(fp_rmin))
        stat.set('input', 'fp_rmax', '{}'.format(fp_rmax))
        stat.set('input', 'fp_npoints', '{}'.format(fp_npoints))
        stat.set('input', 'fp_sigma', '{}'.format(fp_sigma))
        stat.set('input', 'maxgen', '{}'.format(maxgen))

    # ---------- LAQA
    if algo == 'LAQA':
        stat.set('input', 'nselect', '{}'.format(nselect))
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
        stat.set('input', 'maxgen', '{}'.format(maxgen))

    # ---------- lattice
    stat.set('input', 'minlen', '{}'.format(minlen))
    stat.set('input', 'maxlen', '{}'.format(maxlen))
    stat.set('input', 'dangle', '{}'.format(dangle))
    for i in range(len(atype)):
        stat.set('input', 'mindist_{}'.format(i+1), '{}'.format(' '.join(str(c) for c in mindist[i])))

    # ---------- VASP
    if calc_code == 'VASP':
        stat.set('input', 'kppvol', '{}'.format(' '.join(str(c) for c in kppvol)))
        stat.set('input', 'force_gamma', '{}'.format(force_gamma))

    # ---------- QE
    if calc_code == 'QE':
        stat.set('input', 'qe_infile', '{}'.format(qe_infile))
        stat.set('input', 'qe_outfile', '{}'.format(qe_outfile))
        stat.set('input', 'kppvol', '{}'.format(' '.join(str(c) for c in kppvol)))
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
        stat.set('input', 'lammps_potential', '{}'.format(lammps_potential))
        stat.set('input', 'lammps_data', '{}'.format(lammps_data))

    # ---------- option
    stat.set('input', 'maxcnt', '{}'.format(maxcnt))
    stat.set('input', 'stop_chkpt', '{}'.format(stop_chkpt))
    stat.set('input', 'symprec', '{}'.format(symprec))
    if spgnum == 0 or spgnum == 'all':
        stat.set('input', 'spgnum', '{}'.format(spgnum))
    else:
        stat.set('input', 'spgnum', '{}'.format(' '.join(str(d) for d in spgnum)))
    stat.set('input', 'load_struc_flag', '{}'.format(load_struc_flag))
    stat.set('input', 'stop_next_struc', '{}'.format(stop_next_struc))
    stat.set('input', 'append_struc_ea', '{}'.format(append_struc_ea))
    stat.set('input', 'energy_step_flag', '{}'.format(energy_step_flag))
    stat.set('input', 'struc_step_flag', '{}'.format(struc_step_flag))
    stat.set('input', 'fs_step_flag', '{}'.format(fs_step_flag))

    # ---------- write stat
    with open('cryspy.stat', 'w') as fstat:
        stat.write(fstat)


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
    old_nat = [int(x) for x in old_nat.split()]    # character --> integer list
    old_nstage = stat.getint('input', 'nstage')
    old_njob = stat.getint('input', 'njob')
    old_jobcmd = stat.get('input', 'jobcmd')
    old_jobfile = stat.get('input', 'jobfile')

    # ------ BO
    if old_algo == 'BO':
        old_interval = stat.getint('input', 'interval')
        old_score = stat.get('input', 'score')
        old_num_rand_basis = stat.getint('input', 'num_rand_basis')
        old_cdev = stat.getfloat('input', 'cdev')
        old_dscrpt = stat.get('input', 'dscrpt')
        old_fp_rmin = stat.getfloat('input', 'fp_rmin')
        old_fp_rmax = stat.getfloat('input', 'fp_rmax')
        old_fp_npoints = stat.getint('input', 'fp_npoints')
        old_fp_sigma = stat.getfloat('input', 'fp_sigma')
        old_maxgen = stat.getint('input', 'maxgen')

    # ------ LAQA
    if old_algo == 'LAQA':
        old_nselect = stat.getint('input', 'nselect')
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
        old_maxgen = stat.getint('input', 'maxgen')
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
        old_kppvol = [int(x) for x in old_kppvol.split()]    # character --> int
        old_force_gamma = stat.getboolean('input', 'force_gamma')

    # ------ QE
    if old_calc_code == 'QE':
        old_qe_infile = stat.get('input', 'qe_infile')
        old_qe_outfile = stat.get('input', 'qe_outfile')
        old_kppvol = stat.get('input', 'kppvol')
        old_kppvol = [int(x) for x in old_kppvol.split()]    # character --> int
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
        old_spgnum = [int(x) for x in old_spgnum.split()]    # character --> integer list
    old_load_struc_flag = stat.getboolean('input', 'load_struc_flag')
    old_stop_next_struc = stat.getboolean('input', 'stop_next_struc')
    old_append_struc_ea = stat.getboolean('input', 'append_struc_ea')
    old_energy_step_flag = stat.getboolean('input', 'energy_step_flag')
    old_struc_step_flag = stat.getboolean('input', 'struc_step_flag')
    old_fs_step_flag = stat.getboolean('input', 'fs_step_flag')

    # ---------- check difference
    # ------ basic
    if not old_algo == algo:
        raise ValueError('Do not change algo')
        logic_change = True
    if not old_calc_code == calc_code:
        raise ValueError('Do not change calc code')
        logic_change = True
    if not old_tot_struc == tot_struc:
        if algo == 'EA':
            raise ValueError('Do not change tot_struc in EA')
        print('Changed tot_struc from {0} to {1}'.format(old_tot_struc, tot_struc))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed tot_struc from {0} to {1}\n'.format(old_tot_struc, tot_struc))
        logic_change = True
    if not old_natot == natot:
        raise ValueError('Do not change natot')
    if not old_atype == atype:
        raise ValueError('Do not change atype')
    if not old_nat == nat:
        raise ValueError('Do not change nat')
    if not old_nstage == nstage:
        print('Changed nstage from {0} to {1}'.format(old_nstage, nstage))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed nstage from {0} to {1}\n'.format(old_nstage, nstage))
        logic_change = True
    if not old_njob == njob:
        print('Changed njob from {0} to {1}'.format(old_njob, njob))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed njob from {0} to {1}\n'.format(old_njob, njob))
        logic_change = True
    if not old_jobcmd == jobcmd:
        print('Changed jobcmd from {0} to {1}'.format(old_jobcmd, jobcmd))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed jobcmd from {0} to {1}\n'.format(old_jobcmd, jobcmd))
        logic_change = True
    if not old_jobfile == jobfile:
        print('Changed jobfile from {0} to {1}'.format(old_jobfile, jobfile))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed jobfile from {0} to {1}\n'.format(old_jobfile, jobfile))
        logic_change = True

    # ------ BO
    if algo == 'BO':
        if not old_interval == interval:
            print('Changed interval from {0} to {1}'.format(old_interval, interval))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed interval from {0} to {1}\n'.format(old_interval, interval))
                fout.write('####     This will be enabled in next generation\n')
            logic_change = True
        if not old_score == score:
            print('Changed score from {0} to {1}'.format(old_score, score))
            print('    This will be enabled in next BO')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed score from {0} to {1}\n'.format(old_score, score))
                fout.write('####     This will be enabled in next BO\n')
            logic_change = True
        if not old_num_rand_basis == num_rand_basis:
            print('Changed num_rand_basis from {0} to {1}'.format(old_num_rand_basis, num_rand_basis))
            print('    This will be enabled in next BO')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed num_rand_basis from {0} to {1}\n'.format(old_num_rand_basis, num_rand_basis))
                fout.write('####     This will be enabled in next BO\n')
            logic_change = True
        if not old_cdev == cdev:
            print('Changed cdev from {0} to {1}'.format(old_cdev, cdev))
            print('    This will be enabled in next BO')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed cdev from {0} to {1}\n'.format(old_cdev, cdev))
                fout.write('####     This will be enabled in next BO\n')
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
        if not old_maxgen == maxgen:
            print('Changed maxgen from {0} to {1}'.format(old_maxgen, maxgen))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed maxgen from {0} to {1}\n'.format(old_maxgen, maxgen))
            logic_change = True

    # ------ LAQA
    if algo == 'LAQA':
        if not old_nselect == nselect:
            print('Changed nselect from {0} to {1}'.format(old_nselect, nselect))
            print('    This will be enabled in next selection')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed nselect from {0} to {1}\n'.format(old_nselect, nselect))
                fout.write('####     This will be enabled in next selection\n')
            logic_change = True
        if not old_weight_laqa == weight_laqa:
            print('Changed weight_laqa from {0} to {1}'.format(old_weight_laqa, weight_laqa))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed weight_laqa from {0} to {1}\n'.format(old_weight_laqa, weight_laqa))
            logic_change = True

    # ------ EA
    if algo == 'EA':
        if not old_n_pop == n_pop:
            print('Changed n_pop from {0} to {1}'.format(old_n_pop, n_pop))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_pop from {0} to {1}\n'.format(old_n_pop, n_pop))
            logic_change = True
        if not old_n_crsov == n_crsov:
            print('Changed n_crsov from {0} to {1}'.format(old_n_crsov, n_crsov))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_crsov from {0} to {1}\n'.format(old_n_crsov, n_crsov))
            logic_change = True
        if not old_n_perm == n_perm:
            print('Changed n_perm from {0} to {1}'.format(old_n_perm, n_perm))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_perm from {0} to {1}\n'.format(old_n_perm, n_perm))
            logic_change = True
        if not old_n_strain == n_strain:
            print('Changed n_strain from {0} to {1}'.format(old_n_strain, n_strain))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_n_strain from {0} to {1}\n'.format(old_n_strain, n_strain))
            logic_change = True
        if not old_n_rand == n_rand:
            print('Changed n_rand from {0} to {1}'.format(old_n_rand, n_rand))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_fittest from {0} to {1}\n'.format(old_n_rand, n_rand))
            logic_change = True
        if not old_n_elite == n_elite:
            print('Changed n_elite from {0} to {1}'.format(old_n_elite, n_elite))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_elite from {0} to {1}\n'.format(old_n_elite, n_elite))
            logic_change = True
        if not old_fit_reverse == fit_reverse:
            raise ValueError('Do not change fit_reverse')
        if not old_n_fittest == n_fittest:
            print('Changed n_fittest from {0} to {1}'.format(old_n_fittest, n_fittest))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed n_fittest from {0} to {1}\n'.format(old_n_fittest, n_fittest))
            logic_change = True
        if not old_slct_func == slct_func:
            print('Changed slct_func from {0} to {1}'.format(old_slct_func, slct_func))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed slct_func from {0} to {1}\n'.format(old_slct_func, slct_func))
            logic_change = True
        if old_slct_func == 'TNM' and slct_func == 'TNM':
            if not old_t_size == t_size:
                print('Changed t_size from {0} to {1}'.format(old_t_size, t_size))
                print('    This will be enabled in next generation')
                with open('cryspy.out', 'a') as fout:
                    fout.write('\n### Changed t_size from {0} to {1}\n'.format(old_t_size, t_size))
                logic_change = True
        elif old_slct_func == 'RLT' and slct_func == 'RLT':
            if not old_a_rlt == a_rlt:
                print('Changed a_rlt from {0} to {1}'.format(old_a_rlt, a_rlt))
                print('    This will be enabled in next generation')
                with open('cryspy.out', 'a') as fout:
                    fout.write('\n### Changed a_rlt from {0} to {1}\n'.format(old_a_rlt, a_rlt))
                logic_change = True
            if not old_b_rlt == b_rlt:
                print('Changed b_rlt from {0} to {1}'.format(old_b_rlt, b_rlt))
                print('    This will be enabled in next generation')
                with open('cryspy.out', 'a') as fout:
                    fout.write('\n### Changed b_rlt from {0} to {1}\n'.format(old_b_rlt, b_rlt))
                logic_change = True
        elif not old_slct_func == slct_func:
            print('Changed slct_func from {0} to {1}'.format(old_slct_func, slct_func))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed slct_func from {0} to {1}\n'.format(old_slct_func, slct_func))
            logic_change = True
        if not old_crs_lat == crs_lat:
            print('Changed crs_lat from {0} to {1}'.format(old_crs_lat, crs_lat))
            print('    This will be enabled in next generation')
            with open('CrySPY.out', 'a') as fout:
                fout.write('\n### Changed crs_lat from {0} to {1}\n'.format(old_crs_lat, crs_lat))
            logic_change = True
        if not old_crs_func == crs_func:
            print('Changed crs_func from {0} to {1}'.format(old_crs_func, crs_func))
            print('    This will be enabled in next generation')
            with open('cryspy.out', 'a') as fout:
                fout.write('\n### Changed crs_func from {0} to {1}\n'.format(old_crs_func, crs_func))
            logic_change = True
        if not old_nat_diff_tole == nat_diff_tole:
            print('Changed nat_diff_tole from {0} to {1}'.format(old_nat_diff_tole, nat_diff_tole))
            print('    This will be enabled in next generation')
            with open('CrySPY.out', 'a') as fout:
                fout.write('\n### Changed nat_diff_tole from {0} to {1}\n'.format(old_nat_diff_tole, nat_diff_tole))
            logic_change = True
        if not old_ntimes == ntimes:
            print('Changed ntimes from {0} to {1}'.format(old_ntimes, ntimes))
            print('    This will be enabled in next generation')
            with open('CrySPY.out', 'a') as fout:
                fout.write('\n### Changed ntimes from {0} to {1}\n'.format(old_ntimes, ntimes))
            logic_change = True
        if not old_sigma_st == sigma_st:
            print('Changed sigma_st from {0} to {1}'.format(old_sigma_st, sigma_st))
            print('    This will be enabled in next generation')
            with open('CrySPY.out', 'a') as fout:
                fout.write('\n### Changed sigma_st from {0} to {1}\n'.format(old_sigma_st, sigma_st))
            logic_change = True
        if not old_maxcnt_ea == maxcnt_ea:
            print('Changed maxcnt_ea from {0} to {1}'.format(old_maxcnt_ea, maxcnt_ea))
            print('    This will be enabled in next generation')
            with open('CrySPY.out', 'a') as fout:
                fout.write('\n### Changed maxcnt_ea from {0} to {1}\n'.format(old_maxcnt_ea, maxcnt_ea))
            logic_change = True
        if not old_maxgen == maxgen:
            print('Changed maxgen from {0} to {1}'.format(old_maxgen, maxgen))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed maxgen from {0} to {1}\n'.format(old_maxgen, maxgen))
            logic_change = True

    # ------ lattice
    if not old_minlen == minlen:
        print('Changed minlen from {0} to {1}'.format(old_minlen, minlen))
        print('    This will be enabled in next structure generation')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed minlen from {0} to {1}\n'.format(old_minlen, minlen))
            fout.write('####     This will be enabled in next structure generation\n')
        logic_change = True
    if not old_maxlen == maxlen:
        print('Changed maxlen from {0} to {1}'.format(old_maxlen, maxlen))
        print('    This will be enabled in next structure generation')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed maxlen from {0} to {1}\n'.format(old_maxlen, maxlen))
            fout.write('####     This will be enabled in next structure generation\n')
        logic_change = True
    if not old_dangle == dangle:
        print('Changed dangle from {0} to {1}'.format(old_dangle, dangle))
        print('    This will be enabled in next structure generation')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed dangle from {0} to {1}\n'.format(old_dangle, dangle))
            fout.write('####     This will be enabled in next structure generation\n')
        logic_change = True
    if not old_mindist == mindist:
        print('Changed mindist from {0} to {1}'.format(old_mindist, mindist))
        print('    This will be enabled in next structure generation')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed mindist from {0} to {1}\n'.format(old_mindist, mindist))
            fout.write('####     This will be enabled in next structure generation\n')
        logic_change = True

    # ------ VASP
    if calc_code == 'VASP':
        if not old_kppvol == kppvol:
            print('Changed kppvol from {0} to {1}'.format(old_kppvol, kppvol))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed kppvol from {0} to {1}\n'.format(old_kppvol, kppvol))
            logic_change = True
        if not old_force_gamma == force_gamma:
            print('Changed force_gamma from {0} to {1}'.format(old_force_gamma, force_gamma))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed force_gamma from {0} to {1}\n'.format(old_force_gamma, force_gamma))
            logic_change = True

    # ------ QE
    if calc_code == 'QE':
        if not old_qe_infile == qe_infile:
            raise ValueError('Do not change qe_infile')
        if not old_qe_outfile == qe_outfile:
            raise ValueError('Do not change qe_outfile')
        if not old_kppvol == kppvol:
            print('Changed kppvol from {0} to {1}'.format(old_kppvol, kppvol))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed kppvol from {0} to {1}\n'.format(old_kppvol, kppvol))
            logic_change = True
        if not old_force_gamma == force_gamma:
            print('Changed force_gamma from {0} to {1}'.format(old_force_gamma, force_gamma))
            with open('cryspy.out', 'a') as fout:
                fout.write('\n#### Changed force_gamma from {0} to {1}\n'.format(old_force_gamma, force_gamma))
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
        print('Changed maxcnt from {0} to {1}'.format(old_maxcnt, maxcnt))
        print('    This will be enabled in next structure generation')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed maxcnt from {0} to {1}\n'.format(old_maxcnt, maxcnt))
            fout.write('####     This will be enabled in next structure generation\n')
        logic_change = True
    if not old_stop_chkpt == stop_chkpt:
        print('Changed stop_chkpt from {0} to {1}'.format(old_stop_chkpt, stop_chkpt))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed stop_chkpt from {0} to {1}\n'.format(old_stop_chkpt, stop_chkpt))
        logic_change = True
    if not old_symprec == symprec:
        print('Changed symprec from {0} to {1}'.format(old_symprec, symprec))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed symprec from {0} to {1}\n'.format(old_symprec, symprec))
        logic_change = True
    if not old_spgnum == spgnum:
        print('Changed spgnum')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed spgnum}\n')
        logic_change = True
    if not old_load_struc_flag == load_struc_flag:
        print('Changed load_struc_flag from {0} to {1}'.format(old_load_struc_flag, load_struc_flag))
        print('    load_struc_flag is useless when you continue simulations')
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed load_struc_flag from {0} to {1}\n'.format(old_load_struc_flag, load_struc_flag))
            fout.write('####    load_struc_flag is useless when you continue simulations\n')
        logic_change = True
    if not old_stop_next_struc == stop_next_struc:
        print('Changed stop_next_struc from {0} to {1}'.format(old_stop_next_struc, stop_next_struc))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed stop_next_struc from {0} to {1}\n'.format(old_stop_next_struc, stop_next_struc))
        logic_change = True
    if not old_energy_step_flag == energy_step_flag:
        raise ValueError('Do not change energy_step_flag')
    if not old_struc_step_flag == struc_step_flag:
        raise ValueError('Do not change struc_step_flag')
    if not old_fs_step_flag == fs_step_flag:
        raise ValueError('Do not change fs_step_flag')

    # ---------- save stat if necessary
    if logic_change:
        save_stat(stat)
