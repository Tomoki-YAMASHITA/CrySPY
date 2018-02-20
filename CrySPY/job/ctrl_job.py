#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import subprocess

import numpy as np
import pandas as pd

from ..interface import select_code
from ..IO import read_input as rin
from ..IO import pkl_data
from ..IO import out_struc
from ..IO.out_results import out_rslt
from ..IO.out_results import out_LAQA_status, out_LAQA_step, out_LAQA_score
from ..IO.out_results import out_LAQA_energy, out_LAQA_bias, out_LAQA_id_hist
from ..BO import combo_cryspy
from ..BO import select_descriptor
from ..LAQA.calc_score import calc_LAQA_bias


class Ctrl_job(object):

    def __init__(self, stat, init_struc_data, opt_struc_data, rslt_data):
        self.stat = stat
        self.init_struc_data = init_struc_data
        self.opt_struc_data = opt_struc_data
        self.rslt_data = rslt_data

    def RS_init(self, RS_id_data):
        self.next_id, self.id_done = RS_id_data

    def BO_init(self, BO_id_data, BO_data):
        self.gen, self.non_error_id, self.id_to_calc, self.id_done = BO_id_data
        self.descriptors, self.targets = BO_data
        self.logic_next_gen = False

    def LAQA_init(self, LAQA_id_data, LAQA_data):
        self.id_to_calc, self.id_select_hist, self.id_done = LAQA_id_data
        (self.tot_step_select, self.LAQA_step, self.LAQA_struc,
         self.LAQA_energy, self.LAQA_bias, self.LAQA_score) = LAQA_data
        self.logic_next_selection = False

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
        print('work{0:04d}: Structure ID {1} Stage {2} Done!'.format(
            self.work_id, self.id_stat[self.work_id], self.stage_stat[self.work_id]))
        # ---------- next stage
        if self.stage_stat[self.work_id] < rin.nstage:
            self.ctrl_next_stage()
        # ---------- collect result and next struc
        elif self.stage_stat[self.work_id] == rin.nstage:
            self.ctrl_collect()
            self.ctrl_next_struc()
        # ---------- error
        else:
            raise ValueError('Wrong stage in '+self.work_path+'stat_job')

    def ctrl_next_stage(self):
        current_id = self.id_stat[self.work_id]
        # ---------- energy step
        if rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(self.energy_step_data, current_id, self.work_path)
        # ---------- struc step
        if rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(self.struc_step_data, current_id, self.work_path)
        # ---------- fs step
        if rin.fs_step_flag:
            self.fs_step_data = select_code.get_fs_step(self.fs_step_data, current_id, self.work_path)
        # ---------- next stage
        if rin.kpt_flag:
            skip_flag, self.kpt_data = select_code.next_stage(self.stage_stat[self.work_id] + 1,
                                                              self.work_path, self.kpt_data, current_id)
        else:
            skip_flag = select_code.next_stage(self.stage_stat[self.work_id] + 1,
                                               self.work_path)
        # ---------- skip
        if skip_flag:
            self.ctrl_skip()
            return
        # ---------- prepare jobfile
        self.prepare_jobfile(current_id)
        # ---------- submit
        self.submit_next_stage()

    def submit_next_stage(self):
        # ---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<8}    # Structure ID\n'.format(self.id_stat[self.work_id]))
            fwstat.write('{:<8}    # Stage\n'.format(self.stage_stat[self.work_id] + 1))
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([rin.jobcmd, rin.jobfile], stdout=logf, stderr=logf)
        os.chdir('../')    # go back to ..
        # ---------- save status: work???? = structure_ID, Stage
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage {1}'.format(self.id_stat[self.work_id],
                                                    self.stage_stat[self.work_id] + 1))
        with open('cryspy.stat', 'w') as fstat:
            self.stat.write(fstat)

        print('    submit job, structure ID {0} Stage {1}'.format(
            self.id_stat[self.work_id], self.stage_stat[self.work_id] + 1))

    def ctrl_collect(self):
        # ---------- collect results
        current_id = self.id_stat[self.work_id]
        # ---------- energy step
        if rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(self.energy_step_data, current_id, self.work_path)
        # ---------- struc step
        if rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(self.struc_step_data, current_id, self.work_path)
        # ---------- fs step
        if rin.fs_step_flag:
            self.fs_step_data = select_code.get_fs_step(self.fs_step_data, current_id, self.work_path)
        # ---------- each algo
        if rin.algo == 'RS':
            self.ctrl_collect_RS(current_id)
        elif rin.algo == 'BO':
            self.ctrl_collect_BO(current_id)
        elif rin.algo == 'LAQA':
            self.ctrl_collect_LAQA(current_id)
        else:
            raise ValueError('Error, algo')

    def ctrl_collect_RS(self, current_id):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(current_id, self.work_path)
        with open('cryspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(current_id, energy))
        print('    collect results: E = {0}'.format(energy))
        # ---------- get initial spg info
        spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)
        # ---------- success
        if opt_struc is not None:
            # ------ get opt spg info
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symtoleR)
            # ------ out opt_struc
            out_struc.out_opt_struc(opt_struc, current_id)
            out_struc.out_opt_cif(opt_struc, current_id, self.work_path)
        # ---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = None
        # ---------- register opt_struc
        self.opt_struc_data[current_id] = opt_struc
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- save rslt
        tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                energy, magmom, check_opt], index=self.rslt_data.columns)
        self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ------ success
        if opt_struc is not None:
            # -- register id_done
            self.id_done = np.r_[self.id_done, np.array([current_id])]
            # -- save
            RS_id_data = (self.next_id, self.id_done)
            pkl_data.save_RS_id(RS_id_data)

    def ctrl_collect_BO(self, current_id):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(current_id, self.work_path)
        with open('cryspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(current_id, energy))
        print('    collect results: E = {0}'.format(energy))
        # ---------- get initial spg info
        spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)
        # ---------- success
        if opt_struc is not None:
            # ------ get opt spg info
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symtoleR)
            # ------ out opt_struc
            out_struc.out_opt_struc(opt_struc, current_id)
            out_struc.out_opt_cif(opt_struc, current_id, self.work_path)
        # ---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = None
        # ---------- register opt_struc
        self.opt_struc_data[current_id] = opt_struc
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- save rslt
        tmp_series = pd.Series([self.gen, current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                energy, magmom, check_opt], index=self.rslt_data.columns)
        self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ---------- index
        neid_index = np.where(self.non_error_id == current_id)[0]    # np.where returns tuple
        # ---------- success
        if opt_struc is not None:
            # ------ calc descriptor for opt sturcture
            dscrpt = select_descriptor.calc_X([opt_struc])
            # ------ register id_done and targets
            self.id_done = np.r_[self.id_done, np.array([current_id])]
            self.targets = np.r_[self.targets, np.array([energy])]
            # ------ replace opt_descriptor
            self.descriptors[neid_index[0]] = dscrpt     # neid_index[0]: int
        # ---------- error
        else:
            # ------ remove data
            self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
            self.descriptors = np.delete(self.descriptors, neid_index[0], 0)
        # ---------- save
        BO_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
        pkl_data.save_BO_id(BO_id_data)
        BO_data = (self.descriptors, self.targets)
        pkl_data.save_BO_data(BO_data)

    def ctrl_collect_LAQA(self, current_id):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(current_id, self.work_path, check_file='OUTCAR')
        # ---------- total step and LAQA_step
        #     fs_step_data[0] <-- force_step_data
        #     force_step_data[key][stage][step][atom]
        if self.fs_step_data[0][current_id][-1] is None:
            self.LAQA_step[current_id].append(0)
        else:
            self.tot_step_select[-1] += len(self.fs_step_data[0][current_id][-1])
            self.LAQA_step[current_id].append(len(self.fs_step_data[0][current_id][-1]))
        # ------ save status
        self.stat.set('status', 'total step', '{}'.format(sum(self.tot_step_select)))
        with open('cryspy.stat', 'w') as fstat:
                self.stat.write(fstat)
        # ---------- append LAQA struc
        self.LAQA_struc[current_id].append(opt_struc)
        # ---------- append LAQA energy
        self.LAQA_energy[current_id].append(energy/rin.natot)
        # ---------- append LAQA bias
        #     fs_step_data[0] <-- force_step_data
        #     force_step_data[key][stage][step][atom]
        tmp_LAQA_bias = calc_LAQA_bias(self.fs_step_data[0][current_id][-1], c=rin.weight_LAQA)
        self.LAQA_bias[current_id].append(tmp_LAQA_bias)
        # ---------- append LAQA score
        if check_opt is 'done':
            self.LAQA_score[current_id].append(-float('inf'))
        elif np.isnan(energy) or np.isnan(tmp_LAQA_bias):
            self.LAQA_score[current_id].append(-float('inf'))
        else:
            self.LAQA_score[current_id].append(-energy/rin.natot + tmp_LAQA_bias)
        # ---------- save LAQA data
        LAQA_data = (self.tot_step_select, self.LAQA_step, self.LAQA_struc,
                     self.LAQA_energy, self.LAQA_bias, self.LAQA_score)
        pkl_data.save_LAQA_data(LAQA_data)
        # ---------- out LAQA data
        out_LAQA_status(self.LAQA_step, self.LAQA_score, self.LAQA_energy, self.LAQA_bias)
        out_LAQA_step(self.LAQA_step)
        out_LAQA_score(self.LAQA_score)
        out_LAQA_energy(self.LAQA_energy)
        out_LAQA_bias(self.LAQA_bias)
        pkl_data.save_LAQA_data(LAQA_data)
        # ---------- case of 'done'
        if check_opt is 'done':
            with open('cryspy.out', 'a') as fout:
                fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(current_id, energy))
            print('    collect results: E = {0}'.format(energy))
            # ------ get initial spg info
            spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)
            # ------ success
            if opt_struc is not None:
                # -- get opt spg info
                spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symtoleR)
                # -- out opt_struc
                out_struc.out_opt_struc(opt_struc, current_id)
                out_struc.out_opt_cif(opt_struc, current_id, self.work_path)
            # ------ error
            else:
                spg_num_opt = 0
                spg_sym_opt = None
            # ------ register opt_struc
            self.opt_struc_data[current_id] = opt_struc
            pkl_data.save_opt_struc(self.opt_struc_data)
            # ------ save rslt
            tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ------ register id_done
            if opt_struc is not None:
                self.id_done.append(current_id)
                # -- save
                LAQA_id_data = (self.id_to_calc, self.id_select_hist, self.id_done)
                pkl_data.save_LAQA_id(LAQA_id_data)

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
                if self.LAQA_struc[self.next_id]:    # vacant list?
                    next_struc_data = self.LAQA_struc[self.next_id][-1]
                else:
                    next_struc_data = self.init_struc_data[self.next_id]
        else:
            raise ValueError('Error, algo')
        # ---------- common part
        if self.next_id < rin.tot_struc:
            print('work{0:04d}: submit job, structure ID {1} Stage 1'.format(
                  self.work_id, self.next_id))
            # ------ prepare input files for structure optimization
            if rin.kpt_flag:
                self.kpt_data = select_code.next_struc(next_struc_data, self.next_id,
                                                       self.work_path, self.kpt_data)
            else:
                select_code.next_struc(next_struc_data, self.next_id, self.work_path)
            # ------ prepare jobfile
            self.prepare_jobfile(self.next_id)
            # ------ submit
            self.submit_next_struc()
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
        # ---------- update status
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage 1'.format(self.next_id))
        # ---------- RS
        if rin.algo == 'RS':
            self.next_id += 1
            RS_id_data = (self.next_id, self.id_done)
            pkl_data.save_RS_id(RS_id_data)
            self.stat.set('status', 'next_id', '{}'.format(self.next_id))
        # ---------- BO
        elif rin.algo == 'BO':
            BO_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_BO_id(BO_id_data)
            self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        # ---------- LAQA
        elif rin.algo == 'LAQA':
            LAQA_id_data = (self.id_to_calc, self.id_select_hist, self.id_done)
            pkl_data.save_LAQA_id(LAQA_id_data)
            if len(self.id_to_calc) > 30:
                self.stat.set('status', 'id_to_calc', '{} IDs'.format(len(self.id_to_calc)))
            else:
                self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))

    def ctrl_skip(self):
        current_id = self.id_stat[self.work_id]
        # ---------- log and out
        with open('cryspy.out', 'a') as fout:
            fout.write('work{0:04d}: Skip Structure ID {1}\n'.format(self.work_id, current_id))
        print('work{0:04d}: Skip Structure ID{1}'.format(self.work_id, current_id))
        # ---------- get initial spg info
        spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)
        # ---------- 'skip' for rslt
        spg_num_opt = 0
        spg_sym_opt = None
        energy = np.nan
        magmom = np.nan
        check_opt = 'skip'
        # ---------- RS
        if rin.algo == 'RS':
            # ------ save rslt
            tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        # ---------- BO
        elif rin.algo == 'BO':
            # ------ save rslt
            tmp_series = pd.Series([self.gen, current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ------ index
            neid_index = np.where(self.non_error_id == current_id)[0]    # np.where returns tuple
            # ------ remove data
            self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
            self.descriptors = np.delete(self.descriptors, neid_index[0], 0)
            # ------ save
            BO_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_BO_id(BO_id_data)
            BO_data = (self.descriptors, self.targets)
            pkl_data.save_BO_data(BO_data)
        # ---------- LAQA
        elif rin.algo == 'LAQA':
            # ------ save rslt
            tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ---------- LAQA data
            self.LAQA_step[current_id].append(0)
            self.LAQA_struc[current_id].append(None)
            self.LAQA_energy[current_id].append(energy)
            self.LAQA_bias[current_id].append(np.nan)
            self.LAQA_score[current_id].append(-float('inf'))
            # ---------- save LAQA data
            LAQA_data = (self.tot_step_select, self.LAQA_step, self.LAQA_struc,
                         self.LAQA_energy, self.LAQA_bias, self.LAQA_score)
            pkl_data.save_LAQA_data(LAQA_data)
            # ---------- out LAQA data
            out_LAQA_status(self.LAQA_step, self.LAQA_score, self.LAQA_energy, self.LAQA_bias)
            out_LAQA_step(self.LAQA_step)
            out_LAQA_score(self.LAQA_score)
            out_LAQA_energy(self.LAQA_energy)
            out_LAQA_bias(self.LAQA_bias)
        # ---------- clean files
        select_code.clean_calc_files(self.work_path)
        # ---------- next struc
        self.ctrl_next_struc()

    def prepare_jobfile(self, current_id):
        if not os.path.isfile('./calc_in/' + rin.jobfile):
            raise IOError('Could not find ./calc_in' + rin.jobfile)
        with open('./calc_in/' + rin.jobfile, 'r') as f:
            lines = f.readlines()
        lines2 = []
        for line in lines:
            lines2.append(line.replace('CrySPY_ID', str(current_id)))
        with open(self.work_path + rin.jobfile, 'w') as f:
            f.writelines(lines2)

    # ---------- BO
    def ctrl_next_gen(self):
        # ------ out and log
        with open('cryspy.out', 'a') as fout:
            fout.write('# ------ Bayesian optimization\n')
        print('# ------ Bayesian optimization')
        # ------ id_done --> sact
        sact = np.array([], dtype=int)
        for i in self.id_done:
            tindx = np.where(self.non_error_id == i)[0][0]
            sact = np.r_[sact, np.array([tindx])]
        # ------ Bayesian optimization
        actions = combo_cryspy.bayes_opt(sact, self.descriptors, self.targets)
        # ------ actions --> id_to_calc
        for i in actions:
            self.id_to_calc = np.r_[self.id_to_calc, self.non_error_id[i]]
        # ------ gen+1
        self.gen += 1
        # ------ save
        BO_id_data = (self.gen, self.non_error_id, self.id_to_calc, self.id_done)
        pkl_data.save_BO_id(BO_id_data)
        # ------ status
        self.stat.set('status', 'generation', '{}'.format(self.gen))
        self.stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        with open('cryspy.stat', 'w') as fstat:
            self.stat.write(fstat)
        # ------ out and log
        print('# ---------- Generation: {}'.format(self.gen))
        print('selected_id: {}'.format(' '.join(str(a) for a in self.id_to_calc)))
        with open('cryspy.out', 'a') as fout:
            fout.write('# ---------- Generation: {}\n'.format(self.gen))
            fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in self.id_to_calc)))

    # ---------- LAQA
    def ctrl_next_selection(self):
        # ------ LAQA selection
        for k, v in sorted(self.LAQA_score.items(), key=lambda x: -x[1][-1]):
            if v[-1] == -float('inf'):
                break
            else:
                self.id_to_calc.append(k)
                if len(self.id_to_calc) == rin.nselect:
                    break
        # ------ done LAQA
        if len(self.id_to_calc) == 0:
            with open('cryspy.out', 'a') as fout:
                fout.write('\nDone LAQA!\n')
            print('\nDone LAQA!')
            raise SystemExit()
        # ------ append id_select_hist and out
        self.id_select_hist.append(self.id_to_calc)
        out_LAQA_id_hist(self.id_select_hist)
        # ------ tot_step_select for next selection
        self.tot_step_select.append(0)
        # ------ save
        LAQA_id_data = (self.id_to_calc, self.id_select_hist, self.id_done)
        pkl_data.save_LAQA_id(LAQA_id_data)
        LAQA_data = (self.tot_step_select, self.LAQA_step, self.LAQA_struc,
                     self.LAQA_energy, self.LAQA_bias, self.LAQA_score)
        pkl_data.save_LAQA_data(LAQA_data)
        # ------ status
        self.stat.set('status', 'LAQA_selection', '{}'.format(len(self.id_select_hist)))
        if len(self.id_to_calc) > 30:
            self.stat.set('status', 'selected_id', '{} IDs'.format(len(self.id_to_calc)))
            self.stat.set('status', 'id_to_calc', '{} IDs'.format(len(self.id_to_calc)))
        else:
            self.stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
            self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        with open('cryspy.stat', 'w') as fstat:
            self.stat.write(fstat)
        # ------ out and log
        print('\n# ---------- LAQA selection {}'.format(len(self.id_select_hist)))
        with open('cryspy.out', 'a') as fout:
            fout.write('\n# ---------- LAQA selection {}\n'.format(len(self.id_select_hist)))
        if len(self.id_to_calc) > 30:
            print('selected_id: {} IDs'.format(len(self.id_to_calc)))
            with open('cryspy.out', 'a') as fout:
                fout.write('selected_id: {} IDs\n\n'.format(len(self.id_to_calc)))
        else:
            print('selected_id: {}\n'.format(' '.join(str(a) for a in self.id_to_calc)))
            with open('cryspy.out', 'a') as fout:
                fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in self.id_to_calc)))
