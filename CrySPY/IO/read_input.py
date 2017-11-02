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
    if algo not in ['RS', 'BO']:
        raise ValueError('algo should be RS or BO')
    calc_code = config.get('basic', 'calc_code')
    if calc_code not in ['VASP', 'QE', 'soiap']:
        raise ValueError('calc_code should be VASP, QE, or soiap for now')
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
            raise ValueError('score should be TS, EI, or PI, check score')
        try:
            num_rand_basis = config.getint('BO', 'num_rand_basis')
        except:
            num_rand_basis = 0
        try:
            cdev = config.getfloat('BO', 'cdev')
        except:
            cdev = 0.001
        dscrpt = config.get('BO', 'dscrpt')
        if dscrpt == 'FP':
            pass
        else:
            raise ValueError('Now FP only')
        # -- parameters for f-fingerprint (optional)
        try:
            fp_rmin = config.getfloat('BO', 'fp_rmin')
        except:
            fp_rmin = 0.5
        try:
            fp_rmax = config.getfloat('BO', 'fp_rmax')
        except:
            fp_rmax = 5.0
        if fp_rmin < 0.0:
            raise ValueError('fp_rmin < 0, check fp_rmin')
        if fp_rmax < fp_rmin:
            raise ValueError('fp_rmax < fp_rmin, check fp_rmin and fp_rmax')
        try:
            fp_npoints = config.getint('BO', 'fp_npoints')
        except:
            fp_npoints = 50
        if fp_npoints <= 0:
            raise ValueError('fp_npoints <= 0, check fp_npoints')
        try:
            fp_sigma = config.getfloat('BO', 'fp_sigma')
        except:
            fp_sigma = 0.2
        if fp_sigma < 0:
            raise ValueError('fp_sigma < 0, check fp_sigma')

    # ---------- lattice
    # ------ global declaration
    global minlen, maxlen, dangle, mindist
    # ------ read intput variables
    minlen = config.getfloat('lattice', 'minlen')
    maxlen = config.getfloat('lattice', 'maxlen')
    dangle = config.getfloat('lattice', 'dangle')
    if dangle < 0.0:
        raise ValueError('dangle < 0.0, dangle should be positive')
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
        except:
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
        except:
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
    else:
        kpt_flag = False
        force_gamma = False
        raise ValueError('calc_code should be VASP, QE, or soiap for now')

    # ---------- option
    # ------ global declaration
    global maxcnt, stop_chkpt, symtoleI, symtoleR, spgnum
    global load_struc_flag, stop_next_struc
    # ------ read intput variables
    try:
        maxcnt = config.getint('option', 'maxcnt')
    except:
        maxcnt = 200
    try:
        stop_chkpt = config.getint('option', 'stop_chkpt')
    except:
        stop_chkpt = 0
    try:
        symtoleI = config.getfloat('option', 'symtoleI')
    except:
        symtoleI = 0.001
    try:
        symtoleR = config.getfloat('option', 'symtoleR')
    except:
        symtoleR = 0.1
    try:
        spgnum = config.get('option', 'spgnum')
        if spgnum == '0':
            spgnum = 0
        else:
            spgnum = spglist(spgnum)
    except:
        spgnum = 'all'
    try:
        load_struc_flag = config.getboolean('option', 'load_struc_flag')
    except:
        load_struc_flag = False
    try:
        stop_next_struc = config.getboolean('option', 'stop_next_struc')
    except:
        stop_next_struc = False


def spglist(spgnum):
    tmpspg = []
    for c in spgnum.split():
        if '-' in c:
            if not len(c.split('-')) == 2:
                raise ValueError('Wrong input in spgnum. ')
            istart = int(c.split('-')[0])
            iend = int(c.split('-')[1])+1
            if istart < 0 or 230 < istart:
                raise ValueError('spgnum should be 1 -- 230')
            if iend < 0 or 231 < iend:
                raise ValueError('spgnum should be 1 -- 230')
            tmpspg += [i for i in range(istart, iend)]
        else:
            if int(c) < 0 or 230 < int(c):
                raise ValueError('spgnum should be 1 -- 230')
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

        # ------ option
        fout.write('# ------ option section\n')
        fout.write('maxcnt = {}\n'.format(maxcnt))
        fout.write('stop_chkpt = {}\n'.format(stop_chkpt))
        fout.write('symtoleI = {}\n'.format(symtoleI))
        fout.write('symtoleR = {}\n'.format(symtoleR))
        if spgnum == 0 or spgnum == 'all':
            fout.write('spgnum = {}\n'.format(spgnum))
        else:
            fout.write('spgnum = {}\n'.format(' '.join(str(d) for d in spgnum)))
        fout.write('load_struc_flag = {}\n'.format(load_struc_flag))
        fout.write('stop_next_struc = {}\n\n\n'.format(stop_next_struc))


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

    # ---------- option
    stat.set('input', 'maxcnt', '{}'.format(maxcnt))
    stat.set('input', 'stop_chkpt', '{}'.format(stop_chkpt))
    stat.set('input', 'symtoleI', '{}'.format(symtoleI))
    stat.set('input', 'symtoleR', '{}'.format(symtoleR))
    if spgnum == 0 or spgnum == 'all':
        stat.set('input', 'spgnum', '{}'.format(spgnum))
    else:
        stat.set('input', 'spgnum', '{}'.format(' '.join(str(d) for d in spgnum)))
    stat.set('input', 'load_struc_flag', '{}'.format(load_struc_flag))
    stat.set('input', 'stop_next_struc', '{}'.format(stop_next_struc))

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

    # ------ option
    old_maxcnt = stat.getint('input', 'maxcnt')
    old_stop_chkpt = stat.getint('input', 'stop_chkpt')
    old_symtoleI = stat.getfloat('input', 'symtoleI')
    old_symtoleR = stat.getfloat('input', 'symtoleR')
    old_spgnum = stat.get('input', 'spgnum')
    if old_spgnum == '0':
        old_spgnum = 0
    elif not old_spgnum == 'all':
        old_spgnum = [int(x) for x in old_spgnum.split()]    # character --> integer list
    old_load_struc_flag = stat.getboolean('input', 'load_struc_flag')
    old_stop_next_struc = stat.getboolean('input', 'stop_next_struc')

    # ---------- check difference
    # ------ basic
    if not old_algo == algo:
        raise ValueError('Do not change algo')
        logic_change = True
    if not old_calc_code == calc_code:
        raise ValueError('Do not change calc code')
        logic_change = True
    if not old_tot_struc == tot_struc:
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
    if not old_symtoleI == symtoleI:
        print('Changed symtoleI from {0} to {1}'.format(old_symtoleI, symtoleI))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed symtoleI from {0} to {1}\n'.format(old_symtoleI, symtoleI))
        logic_change = True
    if not old_symtoleR == symtoleR:
        print('Changed symtoleR from {0} to {1}'.format(old_symtoleR, symtoleR))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n#### Changed symtoleR from {0} to {1}\n'.format(old_symtoleR, symtoleR))
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

    # ---------- save stat if necessary
    if logic_change:
        save_stat(stat)
