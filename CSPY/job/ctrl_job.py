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
from ..IO import out_results
from ..BO import combo_cspy
from ..BO import select_descriptor


class Ctrl_job(object):


    def __init__(self, stat, init_struc_data, opt_struc_data, rslt_data):
        self.stat = stat
        self.init_struc_data = init_struc_data
        self.opt_struc_data = opt_struc_data
        self.rslt_data = rslt_data


    def RS_init(self, RS_id_data):
        self.next_id, self.id_done = RS_id_data


    def BO_init(self, BO_id_data, BO_data):
        self.gen, self.next_BO_id, self.non_error_id, self.id_to_calc, self.id_done = BO_id_data
        self.descriptors, self.targets = BO_data
        self.logic_next_gen = False

    def check_job(self):
        id_stat = []
        stage_stat = []
        job_stat = []
        for i in range(rin.njob):
            stat_path = 'work{:04d}'.format(i) + '/stat_job'
            try:
                with open(stat_path, 'r') as fstat:
                    istat = fstat.readline()    # id
                    sstat = fstat.readline()    # stage
                    jstat = fstat.readline()    # submitted or done or ...
                    id_stat.append(int(istat.split()[0]))
                    stage_stat.append(int(sstat.split()[0]))
                    if jstat[0:3] == 'sub':
                        job_stat.append('submitted')
                    elif jstat[0:4] == 'done':
                        job_stat.append('done')
                    elif jstat[0:4] == 'skip':
                        job_stat.append('skip')
                    else:
                        job_stat.append('else')
            except:
                id_stat.append('no_file')
                stage_stat.append('no_file')
                job_stat.append('no_file')

        self.id_stat = id_stat
        self.stage_stat = stage_stat
        self.job_stat = job_stat


    def handle_done(self):
        print('work{0:04d}: Structure ID {1} Stage {2} Done!'.format(
            self.work_id, self.id_stat[self.work_id], self.stage_stat[self.work_id]))

        #---------- next stage
        if self.stage_stat[self.work_id] < rin.nstage:
            self.ctrl_next_stage()

        #---------- collect result and next struc
        elif self.stage_stat[self.work_id] == rin.nstage:
            self.ctrl_collect()
            self.ctrl_next_struc()
        else:
            raise ValueError('Wrong stage in '+self.work_path+'stat_job')


    def ctrl_next_stage(self):
        print('    submit job, structure ID {0} Stage {1}'.format(
            self.id_stat[self.work_id], self.stage_stat[self.work_id] + 1))
        select_code.next_stage(self.stage_stat[self.work_id] + 1, self.work_path)
        self.submit_next_stage()


    def submit_next_stage(self):
        #---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<8}    # Structure ID\n'.format(self.id_stat[self.work_id]))
            fwstat.write('{:<8}    # Stage\n'.format(self.stage_stat[self.work_id] + 1))
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([rin.jobcmd, rin.jobfile], stdout=logf, stderr=logf)
        os.chdir('../')    # go back to ..

        #----- save status: work???? = structure_ID, Stage
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage {1}'.format(self.id_stat[self.work_id],
                                                    self.stage_stat[self.work_id] + 1))
        with open('cspy.stat', 'w') as fstat:
            self.stat.write(fstat)


    def ctrl_collect(self):
        #---------- collect results
        current_id = self.id_stat[self.work_id]
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(current_id, self.work_path)
        with open('cspy.out', 'a') as fout:
            fout.write('Done! Structure ID {0:>8}: E = {1}\n'.format(current_id, energy))
        print('    collect results: E = {0}'.format(energy))

        #---------- check error
        fail = energy == np.nan or opt_struc == 'error'

        #---------- get initial spg info
        spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)

        #---------- success
        if not fail:
            #------ get opt spg info
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(symprec=rin.symtoleR)

            #------ register opt_struc
            self.opt_struc_data.append(opt_struc)
            pkl_data.save_opt_struc(self.opt_struc_data)
        #---------- error
        else:
            spg_num_opt = 0
            spg_sym_opt = 'error'

        #---------- RS
        if rin.algo == 'RS':
            #------ save rslt
            tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)

            #------ success
            if not fail:
                #------ register id_done
                self.id_done = np.r_[self.id_done, np.array([current_id])]

                #------ save
                RS_id_data = (self.next_id, self.id_done)
                pkl_data.save_RS_id(RS_id_data)


        #---------- BO
        if rin.algo == 'BO':
            #------ save rslt
            tmp_series = pd.Series([self.gen, current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)

            #------ index
            neid_index = np.where(self.non_error_id == current_id)[0]    # np.where returns tuple

            #------ success
            if not fail:
                #-- calc descriptor for opt sturcture
                dscrpt = select_descriptor.calc_X([opt_struc])
                #-- register id_done and targets
                self.id_done = np.r_[self.id_done, np.array([current_id])]
                self.targets = np.r_[self.targets, np.array([energy])]
                #-- replace opt_descriptor
                self.descriptors[neid_index[0]] = dscrpt     # neid_index[0]: int
            #------ error
            else:
                #-- remove data
                self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
                self.descriptors = np.delete(self.descriptors, neid_index[0], 0)

            #------ save
            BO_id_data = (self.gen, self.next_BO_id, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_BO_id(BO_id_data)
            BO_data = (self.descriptors, self.targets)
            pkl_data.save_BO_data(BO_data)

        #---------- save and out rslt
        pkl_data.save_rslt(self.rslt_data)
        out_results.write_rslt(self.rslt_data)


    def ctrl_next_struc(self):
        #---------- option: stop_next_struc
        if rin.stop_next_struc == 1:    # 0:default, 1:stop_next_struc
            #------ clean status in cspy.stat
            self.stat.set('status', 'work{:04d}'.format(self.work_id), '')
            #------ save status
            with open('cspy.stat', 'w') as fstat:
                self.stat.write(fstat)

            return

        #---------- BO
        if rin.algo == 'BO':
            #------ pick up id to calc
            if len(self.id_to_calc) == 0:
                self.logic_next_gen = True
                self.next_id = rin.tot_struc    # to do nothing
            else:
                self.next_id = self.id_to_calc[0]
                self.id_to_calc = self.id_to_calc[1:]

        #---------- common part
        if self.next_id < rin.tot_struc:
            print('work{0:04d}: submit job, structure ID {1} Stage 1'.format(self.work_id, self.next_id))

            #------ prepare input files for structure optimization
            select_code.next_struc(self.init_struc_data, self.next_id, self.work_path)

            #------ prepare jobfile
            if not os.path.isfile('./calc_in/' + rin.jobfile):
                raise IOError('Could not find ./calc_in' + rin.jobfile)
            with open('./calc_in/' + rin.jobfile, 'r') as f:
                lines = f.readlines()
            lines2 = []
            for line in lines:
                lines2.append(line.replace('cspyID', 'cspy' + str(self.next_id)))
            with open(self.work_path + rin.jobfile, 'w') as f:
                f.writelines(lines2)

            #------ submit
            self.submit_next_struc()

        else:
            #------ clean status in cspy.stat
            self.stat.set('status', 'work{:04d}'.format(self.work_id), '')

        #---------- save status
        with open('cspy.stat', 'w') as fstat:
            self.stat.write(fstat)


    def submit_next_struc(self):
        #---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<8}    # Structure ID\n'.format(self.next_id))
            fwstat.write('{:<8}    # Stage\n'.format(1))
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([rin.jobcmd, rin.jobfile], stdout=logf, stderr=logf)
        os.chdir('../')    # go back to csp root dir

        #---------- update status
        self.stat.set('status', 'work{:04d}'.format(self.work_id),
                      'ID {0:>8}, Stage 1'.format(self.next_id))

        #---------- RS
        if rin.algo == 'RS':
            self.next_id += 1
            RS_id_data = (self.next_id, self.id_done)
            pkl_data.save_RS_id(RS_id_data)
            self.stat.set('status', 'next_id', '{}'.format(self.next_id))
        #---------- BO
        elif rin.algo == 'BO':
            BO_id_data = (self.gen, self.next_BO_id, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_BO_id(BO_id_data)
            self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))


    def ctrl_skip(self):
        current_id = self.id_stat[self.work_id]

        #---------- log and out
        with open('cspy.out', 'a') as fout:
            fout.write('work{0:04d}: Skip Structure ID {1}\n'.format(self.work_id, current_id))
        print('work{0:04d}: Skip Structure ID{1}'.format(self.work_id, current_id))

        #---------- get initial spg info
        spg_sym, spg_num = self.init_struc_data[current_id].get_space_group_info(symprec=rin.symtoleI)

        #---------- 'skip' for rslt
        spg_num_opt = 0
        spg_sym_opt = 'skip'
        energy = np.nan
        magmom = np.nan
        check_opt = 'skip'

        #---------- RS
        if rin.algo == 'RS':
            #------ save rslt
            tmp_series = pd.Series([current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_results.write_rslt(self.rslt_data)

        #---------- BO
        elif rin.algo == 'BO':
            #------ save rslt
            tmp_series = pd.Series([self.gen, current_id, spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                                    energy, magmom, check_opt], index=self.rslt_data.columns)
            self.rslt_data = self.rslt_data.append(tmp_series, ignore_index=True)
            pkl_data.save_rslt(self.rslt_data)
            out_results.write_rslt_BO(self.rslt_data)

            #------ index
            neid_index = np.where(self.non_error_id == current_id)[0]    # np.where returns tuple

            #------ remove data
            self.non_error_id = np.delete(self.non_error_id, neid_index[0], 0)
            self.descriptors = np.delete(self.descriptors, neid_index[0], 0)

            #------ save
            BO_id_data = (self.gen, self.next_BO_id, self.non_error_id, self.id_to_calc, self.id_done)
            pkl_data.save_BO_id(BO_id_data)
            BO_data = (self.descriptors, self.targets)
            pkl_data.save_BO_data(BO_data)

        #---------- clean files
        select_code.clean_calc_files(self.work_path)

        #---------- next struc
        self.ctrl_next_struc()


    #---------- BO
    def ctrl_next_gen(self):
        #------ out and log
        with open('cspy.out', 'a') as fout:
            fout.write('\nDone generation {}\n\n'.format(self.gen))
            fout.write('#------ Bayesian optimization\n')
        print('\nDone generation {}\n'.format(self.gen))
        print('#------ Bayesian optimization')

        #------ id_done --> sact
        sact = np.array([], dtype=int)
        for i in self.id_done:
            tindx = np.where(self.non_error_id == i)[0][0]
            sact = np.r_[sact, np.array([tindx])]

        #------ Bayesian optimization
        actions = combo_cspy.bayes_opt(sact, self.descriptors, self.targets)

        #------ actions --> id_to_calc
        for i in actions:
            self.id_to_calc = np.r_[self.id_to_calc, self.non_error_id[i]]

        #------ gen+1
        self.gen += 1

        #------ write and save
        self.BO_save_write()


    #---------- BO
    def BO_save_write(self):
        #------ save
        BO_id_data = (self.gen, self.next_BO_id, self.non_error_id, self.id_to_calc, self.id_done)
        pkl_data.save_BO_id(BO_id_data)

        #------ status
        self.stat.set('status', 'generation', '{}'.format(self.gen))
        self.stat.set('status', 'selected_id', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        self.stat.set('status', 'id_to_calc', '{}'.format(' '.join(str(a) for a in self.id_to_calc)))
        with open('cspy.stat', 'w') as fstat:
            self.stat.write(fstat)

        #------ out and log
        print('#----------Generation: {}'.format(self.gen))
        print('selected_id: {}'.format(' '.join(str(a) for a in self.id_to_calc)))
        with open('cspy.out', 'a') as fout:
            fout.write('#----------Generation: {}\n'.format(self.gen))
            fout.write('selected_id: {}\n\n'.format(' '.join(str(a) for a in self.id_to_calc)))
