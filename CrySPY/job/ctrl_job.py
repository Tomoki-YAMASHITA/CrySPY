#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import subprocess

import numpy as np
import pandas as pd

from ..BO import select_descriptor
from ..gen_struc.struc_util import out_poscar, out_cif
from ..interface import select_code
from ..IO import read_input as rin
from ..IO import pkl_data
from ..IO.out_results import out_rslt
from ..IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
from ..IO.out_results import out_laqa_energy, out_laqa_bias
from ..LAQA.calc_score import calc_laqa_bias


class Ctrl_job(object):

    def __init__(self, stat, init_struc_data, opt_struc_data, rslt_data):
        self.stat = stat
        self.init_struc_data = init_struc_data
        self.opt_struc_data = opt_struc_data
        self.rslt_data = rslt_data

    def rs_init(self, rs_id_data):
        self.next_id, self.id_done = rs_id_data

    def bo_init(self, bo_id_data, bo_data):
        self.gen, self.non_error_id, self.id_to_calc, self.id_done = bo_id_data
        self.descriptors, self.targets = bo_data
        self.logic_next_gen = False

    def laqa_init(self, laqa_id_data, laqa_data):
        self.id_to_calc, self.id_select_hist, self.id_done = laqa_id_data
        (self.tot_step_select, self.laqa_step, self.laqa_struc,
         self.laqa_energy, self.laqa_bias, self.laqa_score) = laqa_data
        self.logic_next_selection = False

    def ea_init(self, ea_id_data):
        self.gen, self.next_id, self.id_done = ea_id_data
        self.logic_next_gen = False

    def check_job(self):
        # ---------- initialize
        self.id_stat = []
        self.stage_stat = []
        self.job_stat = []
        # ---------- check job status
        for i in range(rin.njob):
            stat_path = 'work{:04d}'.format(i) + '/stat_job'
            try:
                with open(stat_path, 'r') as fstat:
                    istat = fstat.readline()    # id
                    sstat = fstat.readline()    # stage
                    jstat = fstat.readline()    # submitted or done or ...
                    self.id_stat.append(int(istat.split()[0]))
                    self.stage_stat.append(int(sstat.split()[0]))
                    if jstat[0:3] == 'sub':
                        self.job_stat.append('submitted')
                    elif jstat[0:4] == 'done':
                        self.job_stat.append('done')
                    elif jstat[0:4] == 'skip':
                        self.job_stat.append('skip')
                    else:
                        self.job_stat.append('else')
            except:
                self.id_stat.append('no_file')
                self.stage_stat.append('no_file')
                self.job_stat.append('no_file')

    def handle_done(self):
        # ---------- log
        print('work{0:04d}: Structure ID {1} Stage {2} Done!'.format(
            self.work_id, self.cID, self.cstage))
        # ---------- next stage
        if self.cstage < rin.nstage:
            self.ctrl_next_stage()
        # ---------- collect result and next struc
        elif self.cstage == rin.nstage:
            self.ctrl_collect()
            self.ctrl_next_struc()
        # ---------- error
        else:
            raise ValueError('Wrong stage in '+self.work_path+'stat_job')

    def ctrl_next_stage(self):
        # ---------- energy step
        if rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(self.energy_step_data, self.cID, self.work_path)
        # ---------- struc step
        if rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(self.struc_step_data, self.cID, self.work_path)
        # ---------- fs step
        if rin.fs_step_flag:
            self.fs_step_data = select_code.get_fs_step(self.fs_step_data, self.cID, self.work_path)
        # ---------- next stage
        if rin.kpt_flag:
            skip_flag, self.kpt_data = select_code.next_stage(self.cstage + 1,
                                                              self.work_path, self.kpt_data, self.cID)
        else:
            skip_flag = select_code.next_stage(self.cstage + 1, self.work_path)
        # ---------- skip
        if skip_flag:
            self.ctrl_skip()
            self.ctrl_next_struc()
            return
        # ---------- prepare jobfile
        self.prepare_jobfile()
        # ---------- submit
        self.submit_next_stage()

    def submit_next_stage(self):
        # ---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<8}    # Structure ID\n'.format(self.cID))
            fwstat.write('{:<8}    # Stage\n'.format(self.cstage + 1))
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([rin.jobcmd, rin.jobfile], stdout=logf, stderr=logf)
        os.chdir('../')    # go back to ..
        # ---------- save status: work???? = structure_ID, Stage
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage {1}'.format(self.cID, self.cstage + 1))
        with open('cryspy.stat', 'w') as fstat:
            self.stat.write(fstat)
        # ---------- log
        print('    submitted job, structure ID {0} Stage {1}'.format(
            self.cID, self.cstage + 1))

    def ctrl_collect(self):
        # ---------- energy step
        if rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(self.energy_step_data, self.cID, self.work_path)
        # ---------- struc step
        if rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(self.struc_step_data, self.cID, self.work_path)
        # ---------- fs step
        if rin.fs_step_flag:
            self.fs_step_data = select_code.get_fs_step(self.fs_step_data, self.cID, self.work_path)
        # ---------- each algo
        if rin.algo == 'RS':
            self.ctrl_collect_rs()
        elif rin.algo == 'BO':
            self.ctrl_collect_bo()
        elif rin.algo == 'LAQA':
            self.ctrl_collect_laqa()
        elif rin.algo == 'EA':
            self.ctrl_collect_ea()
        else:
            raise ValueError('Error, algo')

    def ctrl_collect_rs(self):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.cID, self.work_path)
        with open('cryspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(self.cID, energy))
        print('    collect results: E = {0}'.format(energy))
        # ---------- get initial spg info
        try:
            spg_sym, spg_num = self.init_struc_data[self.cID].get_space_group_info(symprec=rin.symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ---------- success
        if opt_struc is not None:
            # ------ get opt spg info
            try:
                spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symprec)
            except TypeError:
                spg_num_opt = 0
                spg_sym_opt = None
            # ------ out opt_struc
            out_poscar(opt_struc, self.cID, './data/opt_POSCARS')
            try:
                out_cif(opt_struc, self.cID, self.work_path, './data/opt_CIFS.cif', rin.symprec)
            except TypeError:
                print('failed to write opt_CIF')
        # ---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = None
        # ---------- register opt_struc
        self.opt_struc_data[self.cID] = opt_struc
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- save rslt
        tmp_series = pd.Series([self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                energy, magmom, check_opt], index=self.rslt_data.columns)
        self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ------ success
        if opt_struc is not None:
            # -- register id_done
            self.id_done = np.r_[self.id_done, np.array([self.cID])]
            # -- save
            rs_id_data = (self.next_id, self.id_done)
            pkl_data.save_rs_id(rs_id_data)

    def ctrl_collect_bo(self):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.cID, self.work_path)
        with open('cryspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(self.cID, energy))
        print('    collect results: E = {0}'.format(energy))
        # ---------- get initial spg info
        try:
            spg_sym, spg_num = self.init_struc_data[self.cID].get_space_group_info(symprec=rin.symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ---------- success
        if opt_struc is not None:
            # ------ get opt spg info
            try:
                spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symprec)
            except TypeError:
                spg_num_opt = 0
                spg_sym_opt = None
            # ------ out opt_struc
            out_poscar(opt_struc, self.cID, './data/opt_POSCARS')
            try:
                out_cif(opt_struc, self.cID, self.work_path, './data/opt_CIFS.cif', rin.symprec)
            except TypeError:
                print('failed to write opt_CIF')
        # ---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = None
        # ---------- register opt_struc
        self.opt_struc_data[self.cID] = opt_struc
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- save rslt
        tmp_series = pd.Series([self.gen, self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                energy, magmom, check_opt], index=self.rslt_data.columns)
        self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ---------- index
        neid_index = np.where(self.non_error_id == self.cID)[0]    # np.where returns tuple
        # ---------- success
        if opt_struc is not None:
            # ------ calc descriptor for opt sturcture
            dscrpt = select_descriptor.calc_X([opt_struc])
            # ------ register id_done and targets
            self.id_done = np.r_[self.id_done, np.array([self.cID])]
            self.targets = np.r_[self.targets, np.array([energy])]
            # ------ replace opt_descriptor
            self.descriptors[neid_index[0]] = dscrpt     # neid_index[0]: int
        # ---------- error
        else:
            # ------ remove data
            self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
            self.descriptors = np.delete(self.descriptors, neid_index[0], 0)
        # ---------- save
        bo_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
        pkl_data.save_bo_id(bo_id_data)
        bo_data = (self.descriptors, self.targets)
        pkl_data.save_bo_data(bo_data)

    def ctrl_collect_laqa(self):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.cID, self.work_path, check_file='OUTCAR')
        # ---------- total step and laqa_step
        #     fs_step_data[0] <-- force_step_data
        #     force_step_data[key][stage][step][atom]
        if self.fs_step_data[0][self.cID][-1] is None:
            self.laqa_step[self.cID].append(0)
        else:
            self.tot_step_select[-1] += len(self.fs_step_data[0][self.cID][-1])
            self.laqa_step[self.cID].append(len(self.fs_step_data[0][self.cID][-1]))
        # ------ save status
        self.stat.set('status', 'total step', '{}'.format(sum(self.tot_step_select)))
        with open('cryspy.stat', 'w') as fstat:
                self.stat.write(fstat)
        # ---------- append laqa struc
        self.laqa_struc[self.cID].append(opt_struc)
        # ---------- append laqa energy
        self.laqa_energy[self.cID].append(energy/rin.natot)
        # ---------- append laqa bias
        #     fs_step_data[0] <-- force_step_data
        #     force_step_data[key][stage][step][atom]
        tmp_laqa_bias = calc_laqa_bias(self.fs_step_data[0][self.cID][-1], c=rin.weight_laqa)
        self.laqa_bias[self.cID].append(tmp_laqa_bias)
        # ---------- append laqa score
        if check_opt is 'done':
            self.laqa_score[self.cID].append(-float('inf'))
        elif np.isnan(energy) or np.isnan(tmp_laqa_bias):
            self.laqa_score[self.cID].append(-float('inf'))
        else:
            self.laqa_score[self.cID].append(-energy/rin.natot + tmp_laqa_bias)
        # ---------- save laqa data
        laqa_data = (self.tot_step_select, self.laqa_step, self.laqa_struc,
                     self.laqa_energy, self.laqa_bias, self.laqa_score)
        pkl_data.save_laqa_data(laqa_data)
        # ---------- out laqa data
        out_laqa_status(self.laqa_step, self.laqa_score, self.laqa_energy, self.laqa_bias)
        out_laqa_step(self.laqa_step)
        out_laqa_score(self.laqa_score)
        out_laqa_energy(self.laqa_energy)
        out_laqa_bias(self.laqa_bias)
        pkl_data.save_laqa_data(laqa_data)
        # ---------- case of 'done'
        if check_opt is 'done':
            with open('cryspy.out', 'a') as fout:
                fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(self.cID, energy))
            print('    collect results: E = {0}'.format(energy))
            # ------ get initial spg info
            try:
                spg_sym, spg_num = self.init_struc_data[self.cID].get_space_group_info(symprec=rin.symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            # ------ success
            if opt_struc is not None:
                # -- get opt spg info
                try:
                    spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symprec)
                except TypeError:
                    spg_num_opt = 0
                    spg_sym_opt = None
                # -- out opt_struc
                out_poscar(opt_struc, self.cID, './data/opt_POSCARS')
                try:
                    out_cif(opt_struc, self.cID, self.work_path, './data/opt_CIFS.cif', rin.symprec)
                except TypeError:
                    print('failed to write opt_CIF')
            # ------ error
            else:
                spg_num_opt = 0
                spg_sym_opt = None
            # ------ register opt_struc
            self.opt_struc_data[self.cID] = opt_struc
            pkl_data.save_opt_struc(self.opt_struc_data)
            # ------ save rslt
            tmp_series = pd.Series([self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ------ register id_done
            if opt_struc is not None:
                self.id_done.append(self.cID)
                # -- save
                laqa_id_data = (self.id_to_calc, self.id_select_hist, self.id_done)
                pkl_data.save_laqa_id(laqa_id_data)

    def ctrl_collect_ea(self):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.cID, self.work_path)
        with open('cryspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(self.cID, energy))
        print('    collect results: E = {0}'.format(energy))
        # ---------- get initial spg info
        try:
            spg_sym, spg_num = self.init_struc_data[self.cID].get_space_group_info(symprec=rin.symprec)
        except TypeError:
            spg_num = 0
            spg_sym = None
        # ---------- success
        if opt_struc is not None:
            # ------ get opt spg info
            try:
                spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symprec)
            except TypeError:
                spg_num_opt = 0
                spg_sym_opt = None
            # ------ out opt_struc
            out_poscar(opt_struc, self.cID, './data/opt_POSCARS')
            try:
                out_cif(opt_struc, self.cID, self.work_path, './data/opt_CIFS.cif', rin.symprec)
            except TypeError:
                print('failed to write opt_CIF')
        # ---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = None
        # ---------- register opt_struc
        self.opt_struc_data[self.cID] = opt_struc
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- save rslt
        tmp_series = pd.Series([self.gen, self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                energy, magmom, check_opt], index=self.rslt_data.columns)
        self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ------ success
        if opt_struc is not None:
            # -- register id_done
            self.id_done = np.r_[self.id_done, np.array([self.cID])]
            # -- save
            ea_id_data = (self.gen, self.next_id, self.id_done)
            pkl_data.save_ea_id(ea_id_data)

    def ctrl_next_struc(self):
        # ---------- option: stop_next_struc
        if rin.stop_next_struc:
            # ------ clean status in cryspy.stat
            self.stat.set('status', 'work{:04d}'.format(self.work_id), '')
            # ------ save status
            with open('cryspy.stat', 'w') as fstat:
                self.stat.write(fstat)
            # ------ return
            return
        # ---------- RS
        if rin.algo == 'RS':
            if self.next_id < rin.tot_struc:
                next_struc_data = self.init_struc_data[self.next_id]
        # ---------- BO
        elif rin.algo == 'BO':
            # ------ pick up id to calc
            if len(self.id_to_calc) == 0:
                self.logic_next_gen = True
                self.next_id = rin.tot_struc    # to do nothing
            else:
                self.next_id = self.id_to_calc[0]
                self.id_to_calc = self.id_to_calc[1:]
                next_struc_data = self.init_struc_data[self.next_id]
        # ---------- LAQA
        elif rin.algo == 'LAQA':
            # ------ pick up id to calc
            if len(self.id_to_calc) == 0:
                self.logic_next_selection = True
                self.next_id = rin.tot_struc    # to do nothing
            else:
                self.next_id = self.id_to_calc[0]
                self.id_to_calc = self.id_to_calc[1:]
                if self.laqa_struc[self.next_id]:    # vacant list?
                    next_struc_data = self.laqa_struc[self.next_id][-1]
                else:
                    next_struc_data = self.init_struc_data[self.next_id]
        # ---------- EA
        elif rin.algo == 'EA':
            if self.next_id == rin.tot_struc:
                self.logic_next_gen = True
            elif self.next_id < rin.tot_struc:
                next_struc_data = self.init_struc_data[self.next_id]
        # ---------- algo is wrong
        else:
            raise ValueError('Error, algo')
        # ---------- common part
        self.cID = self.next_id
        if self.next_id < rin.tot_struc:
            # ------ in case there is no initial strucure data
            if next_struc_data is None:
                print('work{0:04d}: structure ID {1} is None'.format(
                      self.work_id, self.next_id))
                self.ctrl_skip()
                self.update_status()
            # ------ right initial structure data
            else:
                # -- prepare input files for structure optimization
                if rin.kpt_flag:
                    self.kpt_data = select_code.next_struc(next_struc_data, self.next_id,
                                                           self.work_path, self.kpt_data)
                else:
                    select_code.next_struc(next_struc_data, self.next_id, self.work_path)
                # -- prepare jobfile
                self.prepare_jobfile()
                # -- submit
                print('work{0:04d}: submit job, structure ID {1} Stage 1'.format(
                      self.work_id, self.next_id))
                self.submit_next_struc()
                # -- update status
                self.update_status()
        else:
            # ------ clean status in cryspy.stat
            self.stat.set('status', 'work{:04d}'.format(self.work_id), '')
        # ---------- save status
        with open('cryspy.stat', 'w') as fstat:
            self.stat.write(fstat)

    def submit_next_struc(self):
        # ---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<8}    # Structure ID\n'.format(self.next_id))
            fwstat.write('{:<8}    # Stage\n'.format(1))
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([rin.jobcmd, rin.jobfile], stdout=logf, stderr=logf)
        os.chdir('../')    # go back to csp root dir

    def update_status(self):
        # ---------- update status
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage 1'.format(self.next_id))
        # ---------- RS
        if rin.algo == 'RS':
            self.next_id += 1
            rs_id_data = (self.next_id, self.id_done)
            pkl_data.save_rs_id(rs_id_data)
            self.stat.set('status', 'next_id', '{}'.format(self.next_id))
        # ---------- BO
        elif rin.algo == 'BO':
            bo_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_bo_id(bo_id_data)
            self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        # ---------- LAQA
        elif rin.algo == 'LAQA':
            laqa_id_data = (self.id_to_calc, self.id_select_hist, self.id_done)
            pkl_data.save_laqa_id(laqa_id_data)
            if len(self.id_to_calc) > 30:
                self.stat.set('status', 'id_to_calc', '{} IDs'.format(len(self.id_to_calc)))
            else:
                self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        # ---------- EA
        if rin.algo == 'EA':
            self.next_id += 1
            ea_id_data = (self.gen, self.next_id, self.id_done)
            pkl_data.save_ea_id(ea_id_data)
            self.stat.set('status', 'next_id', '{}'.format(self.next_id))

    def ctrl_skip(self):
        # ---------- log and out
        with open('cryspy.out', 'a') as fout:
            fout.write('work{0:04d}: Skip Structure ID {1}\n'.format(self.work_id, self.cID))
        print('work{0:04d}: Skip Structure ID {1}'.format(self.work_id, self.cID))
        # ---------- get initial spg info
        if self.init_struc_data[self.cID] is None:
            spg_sym = None
            spg_num = 0
        else:
            try:
                spg_sym, spg_num = self.init_struc_data[self.cID].get_space_group_info(symprec=rin.symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
        # ---------- 'skip' for rslt
        spg_num_opt = 0
        spg_sym_opt = None
        energy = np.nan
        magmom = np.nan
        check_opt = 'skip'
        # ---------- register opt_struc
        self.opt_struc_data[self.cID] = None
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- RS
        if rin.algo == 'RS':
            # ------ save rslt
            tmp_series = pd.Series([self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        # ---------- BO
        elif rin.algo == 'BO':
            # ------ save rslt
            tmp_series = pd.Series([self.gen, self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ------ index
            neid_index = np.where(self.non_error_id == self.cID)[0]    # np.where returns tuple
            # ------ remove data
            self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
            self.descriptors = np.delete(self.descriptors, neid_index[0], 0)
            # ------ save
            bo_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_bo_id(bo_id_data)
            bo_data = (self.descriptors, self.targets)
            pkl_data.save_bo_data(bo_data)
        # ---------- LAQA
        elif rin.algo == 'LAQA':
            # ------ save rslt
            tmp_series = pd.Series([self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ---------- laqa data
            self.laqa_step[self.cID].append(0)
            self.laqa_struc[self.cID].append(None)
            self.laqa_energy[self.cID].append(energy)
            self.laqa_bias[self.cID].append(np.nan)
            self.laqa_score[self.cID].append(-float('inf'))
            # ---------- save laqa data
            laqa_data = (self.tot_step_select, self.laqa_step, self.laqa_struc,
                         self.laqa_energy, self.laqa_bias, self.laqa_score)
            pkl_data.save_laqa_data(laqa_data)
            # ---------- out laqa data
            out_laqa_status(self.laqa_step, self.laqa_score, self.laqa_energy, self.laqa_bias)
            out_laqa_step(self.laqa_step)
            out_laqa_score(self.laqa_score)
            out_laqa_energy(self.laqa_energy)
            out_laqa_bias(self.laqa_bias)
        # ---------- EA
        elif rin.algo == 'EA':
            tmp_series = pd.Series([self.gen, self.cID, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        # ---------- clean files
        select_code.clean_calc_files(self.work_path)

    def prepare_jobfile(self):
        if not os.path.isfile('./calc_in/' + rin.jobfile):
            raise IOError('Could not find ./calc_in' + rin.jobfile)
        with open('./calc_in/' + rin.jobfile, 'r') as f:
            lines = f.readlines()
        lines2 = []
        for line in lines:
            lines2.append(line.replace('CrySPY_ID', str(self.cID)))
        with open(self.work_path + rin.jobfile, 'w') as f:
            f.writelines(lines2)
