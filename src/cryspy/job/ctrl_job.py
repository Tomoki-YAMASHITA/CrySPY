import itertools
from logging import getLogger
import os
import shutil
import subprocess

import numpy as np

from ..interface import select_code
from ..IO import change_input, io_stat, pkl_data
from ..IO.out_results import out_rslt, out_hdist
from ..util.utility import backup_cryspy
from ..util.struc_util import out_poscar, out_cif

# ---------- import later
#from ..BO.select_descriptor import select_descriptor
#from ..BO import bo_next_select
#from ..EA import ea_next_gen
#from ..EA.calc_ef import calc_ef
#from ..EA.calc_hull import calc_convex_hull_2d
#from ..LAQA.calc_score import calc_laqa_bias
#from ..LAQA import laqa_next_selection
#from ..IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
#from ..IO.out_results import out_laqa_energy, out_laqa_bias


logger = getLogger('cryspy')


class Ctrl_job:

    def __init__(self, rin, init_struc_data):
        self.rin = rin
        self.init_struc_data = init_struc_data
        self.opt_struc_data = pkl_data.load_opt_struc()
        self.rslt_data = pkl_data.load_rslt()
        self.id_queueing = pkl_data.load_id_queueing()
        self.id_running = pkl_data.load_id_running()
        self.recheck = False
        self.logic_next = False
        # ---------- for each algorithm
        if self.rin.algo == 'RS':
            pass
        elif self.rin.algo == 'BO':
            self.n_selection = pkl_data.load_n_selection()
            self.id_select_hist = pkl_data.load_id_select_hist()
            self.init_dscrpt_data = pkl_data.load_init_dscrpt_data()
            self.opt_dscrpt_data = pkl_data.load_opt_dscrpt_data()
            self.bo_mean = pkl_data.load_bo_mean()
            self.bo_var = pkl_data.load_bo_var()
            self.bo_score = pkl_data.load_bo_score()
        elif self.rin.algo == 'LAQA':
            self.id_select_hist = pkl_data.load_id_select_hist()
            self.tot_step_select = pkl_data.load_tot_step_select()
            self.laqa_step = pkl_data.load_laqa_step()
            self.laqa_struc = pkl_data.load_laqa_struc()
            self.laqa_energy = pkl_data.load_laqa_energy()
            self.laqa_bias = pkl_data.load_laqa_bias()
            self.laqa_score = pkl_data.load_laqa_score()
        elif self.rin.algo in ['EA', 'EA-vc']:
            self.gen = pkl_data.load_gen()
            if self.rin.struc_mode in ['mol', 'mol_bs']:
                self.struc_mol_id = pkl_data.load_struc_mol_id()
            if self.rin.algo == 'EA-vc':
                self.nat_data = pkl_data.load_nat_data()
            # do not have to load ea_data and hdist_data here.
            # ea_data is used only in ea_next_gen.py
            # hdist_data is loaded each time hull distance is calculated
        # ---------- for option
        if self.rin.kpt_flag:
            self.kpt_data = pkl_data.load_kpt()
        if self.rin.energy_step_flag:
            self.energy_step_data = pkl_data.load_energy_step()
        if self.rin.struc_step_flag:
            self.struc_step_data = pkl_data.load_struc_step()
        if self.rin.force_step_flag:
            self.force_step_data = pkl_data.load_force_step()
        if self.rin.stress_step_flag:
            self.stress_step_data = pkl_data.load_stress_step()
        # ---------- flag for next selection or generation
        if self.rin.algo in ['BO', 'LAQA', 'EA', 'EA-vc']:
            if (self.id_queueing or self.id_running):
                self.go_next_sg = False
            else:
                self.go_next_sg = True

    def check_job(self):
        # ---------- option: recalc
        if self.rin.recalc is not None:
            self.set_recalc()
        # ---------- temporarily append
        self.tmp_running = self.id_running[:]    # shallow copy
        self.tmp_queueing = self.id_queueing[:]
        if not self.rin.stop_next_struc:    # if true --> option: stop_next_struc
            while len(self.tmp_running) < self.rin.njob and self.tmp_queueing:
                self.tmp_running.append(self.tmp_queueing.pop(0))
        # ---------- initialize
        self.stage_stat = {}    # key: Structure ID
        self.job_stat = {}
        # ---------- check job status
        for cid in self.tmp_running[:self.rin.njob]:
            # ------ mkdir
            if not os.path.isdir(f'work/{cid}'):
                os.mkdir(f'work/{cid}')
            # ------ check stat_job file
            stat_path = f'work/{cid}' + '/stat_job'
            try:
                with open(stat_path, 'r') as fstat:
                    istat = fstat.readline()    # id
                    sstat = fstat.readline()    # stage
                    jstat = fstat.readline()    # submitted or done or ...
                self.stage_stat[cid] = int(sstat.split()[0])
                if not cid == int(istat.split()[0]):
                    logger.error(f'ID is wrong in work/{cid}')
                    raise SystemExit(1)
                self.stage_stat[cid] = int(sstat.split()[0])
                if jstat[0:3] == 'sub':
                    self.job_stat[cid] = 'submitted'
                elif jstat[0:4] == 'done':
                    self.job_stat[cid] = 'done'
                elif jstat[0:4] == 'skip':
                    self.job_stat[cid] = 'skip'
                else:
                    self.job_stat[cid] = 'else'
            except IOError:
                self.stage_stat[cid] = 'no_file'
                self.job_stat[cid] = 'no_file'

    def set_recalc(self):
        # ---------- check id
        for tid in self.rin.recalc:
            if tid not in self.opt_struc_data:
                logger.error(f'ID {tid} has not yet been calculated')
                raise SystemExit(1)
        # ---------- append IDs to the head of id_queueing
        self.id_queueing = list(self.rin.recalc) + self.id_queueing
        # ---------- save and update id_queueing
        pkl_data.save_id_queueing(self.id_queueing)
        stat = io_stat.stat_read()
        io_stat.set_id(stat, 'id_queueing', self.id_queueing)
        io_stat.write_stat(stat)
        # ---------- log and out
        logger.info('# -- Recalc')
        logger.info(f'Append {self.rin.recalc} to the head of id_queueing')
        # ---------- clear recalc
        self.rin.recalc = None
        config = change_input.read_config()
        change_input.change_input(config, 'option', 'recalc', '')    # clear
        change_input.write_config(config)
        logger.info('Clear recalc in cryspy.in')

    def handle_job(self):
        logger.info('# ---------- job status')
        for cid in self.tmp_running[:self.rin.njob]:
            # ---------- set work_path and cid
            self.work_path = f'./work/{cid}/'
            self.cid = cid
            # ---------- handle job
            if self.job_stat[cid] == 'submitted':
                logger.info(f'ID {cid:>6}: still queueing or running')
            elif self.job_stat[cid] == 'done':
                self.ctrl_done()
            elif self.job_stat[cid] == 'skip':
                self.ctrl_skip()
            elif self.job_stat[cid] == 'else':
                logger.error(f'Wrong job_stat in {self.work_path}. ')
            elif self.job_stat[cid] == 'no_file':
                self.ctrl_next_struc()
            else:
                logger.error(f'Unexpected error in {self.work_path}stat_job')

    def ctrl_done(self):
        self.cstage = self.stage_stat[self.cid]
        # ---------- log
        logger.info(f'ID {self.cid:>6}: Stage {self.cstage} Done!')
        # ---------- next stage
        if self.cstage < self.rin.nstage:
            self.ctrl_next_stage()
        # ---------- collect result
        elif self.cstage == self.rin.nstage:
            self.ctrl_collect()
        # ---------- error
        else:
            logger.error('Wrong stage in '+self.work_path+'stat_job')
            raise SystemExit(1)

    def ctrl_next_stage(self):
        # ---------- EA-vc
        if self.rin.algo == 'EA-vc':
            nat = self.nat_data[self.cid]
        else:
            nat = self.rin.nat
        # ---------- energy step
        if self.rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(
                self.rin, self.energy_step_data, self.cid, self.work_path, nat)
        # ---------- struc step
        if self.rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(
                self.rin, self.struc_step_data, self.cid, self.work_path, nat)
        # ---------- force step
        if self.rin.force_step_flag:
            self.force_step_data = select_code.get_force_step(
                self.rin, self.force_step_data, self.cid, self.work_path, nat)
        # ---------- stress step
        if self.rin.stress_step_flag:
            self.stress_step_data = select_code.get_stress_step(
                self.rin, self.stress_step_data, self.cid, self.work_path)
        # ---------- next stage
        if self.rin.kpt_flag:
            skip_flag, self.kpt_data = select_code.next_stage(
                self.rin, self.cstage, self.work_path, nat,
                kpt_data=self.kpt_data, cid=self.cid)
        else:
            skip_flag = select_code.next_stage(self.rin, self.cstage,
                                               self.work_path, nat)
        # ---------- skip
        if skip_flag:
            self.ctrl_skip()
            return
        # ---------- prepare jobfile
        prepare_jobfile(self.rin, self.cid, self.work_path)
        # ---------- submit
        self.submit_next_stage()

    def submit_next_stage(self):
        # ---------- submit job
        os.chdir(self.work_path)    # cd work_path
        with open('stat_job', 'w') as fwstat:
            fwstat.write(f'{self.cid:<6}    # Structure ID\n')
            fwstat.write(f'{self.cstage + 1:<6}    # Stage\n')
            fwstat.write('submitted\n')
        with open('sublog', 'w') as logf:
            subprocess.Popen([self.rin.jobcmd, self.rin.jobfile],
                             stdout=logf, stderr=logf)
        os.chdir('../../')    # go back to ..
        # ---------- save status
        stat = io_stat.stat_read()
        io_stat.set_stage(stat, self.cid, self.cstage + 1)
        io_stat.write_stat(stat)
        # ---------- log
        logger.info(f'    submitted job, ID {self.cid:>6} Stage {self.cstage + 1}')

    def ctrl_collect(self):
        # ---------- EA-vc
        if self.rin.algo == 'EA-vc':
            nat = self.nat_data[self.cid]
        else:
            nat = self.rin.nat
        # ---------- energy step
        if self.rin.energy_step_flag:
            self.energy_step_data = select_code.get_energy_step(
                self.rin, self.energy_step_data, self.cid, self.work_path, nat)
        # ---------- struc step
        if self.rin.struc_step_flag:
            self.struc_step_data = select_code.get_struc_step(
                self.rin, self.struc_step_data, self.cid, self.work_path, nat)
        # ---------- force step
        if self.rin.force_step_flag:
            self.force_step_data = select_code.get_force_step(
                self.rin, self.force_step_data, self.cid, self.work_path, nat)
        # ---------- stress step
        if self.rin.stress_step_flag:
            self.stress_step_data = select_code.get_stress_step(
                self.rin, self.stress_step_data, self.cid, self.work_path)
        # ---------- each algo
        if self.rin.algo == 'RS':
            self.ctrl_collect_rs(nat)
        elif self.rin.algo == 'BO':
            self.ctrl_collect_bo(nat)
        elif self.rin.algo == 'LAQA':
            self.ctrl_collect_laqa(nat)
        elif self.rin.algo in ['EA', 'EA-vc']:
            self.ctrl_collect_ea(nat)
        else:
            logger.error('Error, algo')
            raise SystemExit(1)
        # ---------- move to fin
        if self.rin.algo == 'LAQA':
            if self.fin_laqa:
                mv_fin(self.cid)
            else:
                os.rename(self.work_path+'stat_job',
                          self.work_path+'prev_stat_job')
        else:
            mv_fin(self.cid)
        # ---------- update status
        self.id_queueing, self.id_running = update_status(
                                                self.cid,
                                                self.id_queueing,
                                                self.id_running,
                                                operation='fin',
                                            )
        # ---------- recheck
        self.recheck = True

    def ctrl_collect_rs(self, nat):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.rin, self.cid, self.work_path, nat)
        logger.info(f'    collect results: E = {energy} eV/atom')
        # ---------- register opt data
        self.opt_struc_data, self.rslt_data = regist_opt(
            self.rin,
            self.cid,
            self.work_path,
            self.init_struc_data,
            self.opt_struc_data,
            self.rslt_data,
            opt_struc,
            energy,
            magmom,
            check_opt,
            ef=None,
            nat=None,
            n_selection=None,
            gen=None,
        )

    def ctrl_collect_bo(self, nat):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.rin, self.cid, self.work_path, nat)
        logger.info(f'    collect results: E = {energy} eV/atom')
        # ---------- register opt data
        self.opt_struc_data, self.rslt_data = regist_opt(
            self.rin,
            self.cid,
            self.work_path,
            self.init_struc_data,
            self.opt_struc_data,
            self.rslt_data,
            opt_struc,
            energy,
            magmom,
            check_opt,
            ef=None,
            nat=None,
            n_selection=self.n_selection,
            gen=None,
        )
        # ---------- success
        if opt_struc is not None:
            from ..BO.select_descriptor import select_descriptor
            # ------ calc descriptor for opt sturcture
            tmp_dscrpt = select_descriptor(self.rin, {self.cid: opt_struc})
            # ------ update descriptors
            self.opt_dscrpt_data.update(tmp_dscrpt)
        # ---------- error
        else:
            # ------ update descriptors and non_error_id
            self.opt_dscrpt_data[self.cid] = None
        # ---------- save bo_data
        pkl_data.save_opt_dscrpt_data(self.opt_dscrpt_data)

    def ctrl_collect_laqa(self, nat):
        # ---------- flag for finish
        self.fin_laqa = False
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.rin, self.cid, self.work_path, nat)
        # ---------- total step and laqa_step
        #     force_step_data[ID][stage][step][atom]
        if self.force_step_data[self.cid][-1] is None:
            self.laqa_step[self.cid].append(0)
        else:
            self.tot_step_select[-1] += len(
                self.force_step_data[self.cid][-1])
            self.laqa_step[self.cid].append(
                len(self.force_step_data[self.cid][-1]))
        # ------ save status
        stat = io_stat.stat_read()
        io_stat.set_common(stat, 'total_step', sum(self.tot_step_select))
        io_stat.write_stat(stat)
        # ---------- append laqa struc
        self.laqa_struc[self.cid].append(opt_struc)
        # ---------- append laqa energy
        self.laqa_energy[self.cid].append(energy)
        # ---------- append laqa bias
        #     force_step_data[ID][stage][step][atom]
        #     stress_step_data[ID][stage][step][atom]
        from ..LAQA.calc_score import calc_laqa_bias
        tmp_laqa_bias = calc_laqa_bias(
            self.force_step_data[self.cid][-1],
            self.stress_step_data[self.cid][-1],
            wf=self.rin.wf, ws=self.rin.ws)
        self.laqa_bias[self.cid].append(tmp_laqa_bias)
        # ---------- append laqa score
        if check_opt == 'done' or np.isnan(energy) or np.isnan(tmp_laqa_bias):
            self.laqa_score[self.cid].append(-float('inf'))
        else:
            self.laqa_score[self.cid].append(-energy + tmp_laqa_bias)
        # ---------- save and out laqa data
        pkl_data.save_tot_step_select(self.tot_step_select)
        pkl_data.save_laqa_step(self.laqa_step)
        pkl_data.save_laqa_struc(self.laqa_struc)
        pkl_data.save_laqa_energy(self.laqa_energy)
        pkl_data.save_laqa_bias(self.laqa_bias)
        pkl_data.save_laqa_score(self.laqa_score)
        from ..IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
        from ..IO.out_results import out_laqa_energy, out_laqa_bias
        out_laqa_status(self.laqa_step, self.laqa_score,
                        self.laqa_energy, self.laqa_bias)
        out_laqa_step(self.laqa_step)
        out_laqa_score(self.laqa_score)
        out_laqa_energy(self.laqa_energy)
        out_laqa_bias(self.laqa_bias)
        # ---------- case of 'done' or error
        if check_opt == 'done' or np.isnan(energy) or np.isnan(tmp_laqa_bias):
            self.fin_laqa = True
            logger.info(f'    collect results: E = {energy} eV/atom')
            # ------ register opt data
            self.opt_struc_data, self.rslt_data = regist_opt(
                self.rin,
                self.cid,
                self.work_path,
                self.init_struc_data,
                self.opt_struc_data,
                self.rslt_data,
                opt_struc,
                energy,
                magmom,
                check_opt,
                ef=None,
                nat=None,
                n_selection=None,
                gen=None,
            )

    def ctrl_collect_ea(self, nat):
        # ---------- get opt data
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(self.rin, self.cid, self.work_path, nat)
        logger.info(f'    collect results: E = {energy} eV/atom')
        # ---------- calculate Ef
        if self.rin.algo == 'EA-vc':
            from ..EA.calc_ef import calc_ef
            ef = calc_ef(energy, nat, self.rin.end_point)
            logger.info(f'                     Ef = {ef} eV/atom')
            regist_nat = nat
        else:
            ef = None
            regist_nat = None
        # ---------- register opt data
        self.opt_struc_data, self.rslt_data = regist_opt(
            self.rin,
            self.cid,
            self.work_path,
            self.init_struc_data,
            self.opt_struc_data,
            self.rslt_data,
            opt_struc,
            energy,
            magmom,
            check_opt,
            ef=ef,
            nat=regist_nat,
            n_selection=None,
            gen=self.gen,
        )

    def ctrl_next_struc(self):
        # ---------- RS, BO, EA, EA-vc
        if self.rin.algo in ['RS', 'BO', 'EA', 'EA-vc']:
            next_struc_data = self.init_struc_data[self.cid]
        # ---------- LAQA
        elif self.rin.algo == 'LAQA':
            if self.laqa_struc[self.cid]:    # vacant list?
                next_struc_data = self.laqa_struc[self.cid][-1]
            else:
                next_struc_data = self.init_struc_data[self.cid]
        # ---------- algo is wrong
        else:
            logger.error('Error, algo in ctrl_next_struc')
            raise SystemExit(1)
        # ---------- nat for EA-vc
        if self.rin.algo == 'EA-vc':
            nat = self.nat_data[self.cid]
        else:
            nat = self.rin.nat
        # ---------- common part
        # ------ in case there is no initial strucure data
        if next_struc_data is None:
            logger.info(f'ID {self.cid:>6}: initial structure is None')
            self.ctrl_skip()
        # ------ normal initial structure data
        else:
            # -- prepare input files for structure optimization
            if self.rin.kpt_flag:
                self.kpt_data = select_code.next_struc(
                    self.rin,
                    next_struc_data,
                    self.cid,
                    self.work_path,
                    nat,
                    self.kpt_data,
                )
            else:
                select_code.next_struc(
                    self.rin,
                    next_struc_data,
                    self.cid,
                    self.work_path,
                    nat,
                )
            # -- prepare jobfile
            prepare_jobfile(self.rin, self.cid, self.work_path)
            # -- submit
            submit_next_struc(self.rin, self.cid, self.work_path)
            logger.info(f'ID {self.cid:>6}: submit job, Stage 1')
            # -- update status
            self.id_queueing, self.id_running = update_status(
                self.cid,
                self.id_queueing,
                self.id_running,
                operation='submit'
            )

    def ctrl_skip(self):
        # ---------- log
        logger.info(f'ID {self.cid:>6}: Skip')
        # ---------- get initial spg info
        if self.init_struc_data[self.cid] is None:
            spg_sym = None
            spg_num = 0
        else:
            try:
                spg_sym, spg_num = self.init_struc_data[
                    self.cid].get_space_group_info(symprec=self.rin.symprec)
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
        self.opt_struc_data[self.cid] = None
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- RS
        if self.rin.algo == 'RS':
            # ------ save rslt
            self.rslt_data.loc[self.cid] = [spg_num, spg_sym,
                                            spg_num_opt, spg_sym_opt,
                                            energy, magmom, check_opt]
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        # ---------- BO
        elif self.rin.algo == 'BO':
            # ------ save rslt
            self.rslt_data.loc[self.cid] = [self.n_selection,
                                            spg_num, spg_sym,
                                            spg_num_opt, spg_sym_opt,
                                            energy, magmom, check_opt]
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ------ update and save descriptors
            self.opt_dscrpt_data[self.cid] = None
            pkl_data.save_opt_dscrpt_data(self.opt_dscrpt_data)
        # ---------- LAQA
        elif self.rin.algo == 'LAQA':
            # ------ save rslt
            self.rslt_data.loc[self.cid] = [spg_num, spg_sym,
                                            spg_num_opt, spg_sym_opt,
                                            energy, magmom, check_opt]
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
            # ---------- laqa data
            self.laqa_step[self.cid].append(0)
            self.laqa_struc[self.cid].append(None)
            self.laqa_energy[self.cid].append(energy)
            self.laqa_bias[self.cid].append(np.nan)
            self.laqa_score[self.cid].append(-float('inf'))
            # ---------- save and out laqa data
            #pkl_data.save_tot_step_select(self.tot_step_select)    # not used here
            pkl_data.save_laqa_step(self.laqa_step)
            pkl_data.save_laqa_struc(self.laqa_struc)
            pkl_data.save_laqa_energy(self.laqa_energy)
            pkl_data.save_laqa_bias(self.laqa_bias)
            pkl_data.save_laqa_score(self.laqa_score)
            from ..IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
            from ..IO.out_results import out_laqa_energy, out_laqa_bias
            out_laqa_status(self.laqa_step, self.laqa_score,
                            self.laqa_energy, self.laqa_bias)
            out_laqa_step(self.laqa_step)
            out_laqa_score(self.laqa_score)
            out_laqa_energy(self.laqa_energy)
            out_laqa_bias(self.laqa_bias)
        # ---------- EA
        elif self.rin.algo == 'EA':
            self.rslt_data.loc[self.cid] = [self.gen,
                                            spg_num, spg_sym,
                                            spg_num_opt, spg_sym_opt,
                                            energy, magmom, check_opt]
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        elif self.rin.algo == 'EA-vc':
            ef = np.nan
            nat = self.nat_data[self.cid]
            self.rslt_data.loc[self.cid] = [self.gen,
                                            spg_num, spg_sym,
                                            spg_num_opt, spg_sym_opt,
                                            energy, ef, nat, magmom, check_opt]
            pkl_data.save_rslt(self.rslt_data)
            out_rslt(self.rslt_data)
        # ---------- move to fin
        mv_fin(self.cid)
        # ---------- update status
        self.id_queueing, self.id_running = update_status(
                                                self.cid,
                                                self.id_queueing,
                                                self.id_running,
                                                operation='fin',
                                            )
        # ---------- recheck
        self.recheck = True

    def next_sg(self, noprint=False):
        '''
        next selection or generation
        '''
        if self.rin.algo == 'BO':
            self.next_select_BO(noprint)
        if self.rin.algo == 'LAQA':
            self.next_select_LAQA()
        if self.rin.algo in ['EA', 'EA-vc']:
            if self.rin.algo == 'EA-vc':
                nat_data = self.nat_data
            else:
                nat_data=None
            next_gen_EA(
                self.rin,
                self.gen,
                self.go_next_sg,
                self.init_struc_data,
                self.opt_struc_data,
                self.rslt_data,
                nat_data=nat_data,
                struc_mol_id=None,
            )

    def next_select_BO(self, noprint=False):
        # ---------- log
        logger.info(f'\nDone selection {self.n_selection}')
        # ---------- done all structures
        if len(self.rslt_data) == self.rin.tot_struc:
            logger.info('\nDone all structures!')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- flag for next selection or generation
        if not self.go_next_sg:
            logger.info('\nBO is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- check point 3
        if self.rin.stop_chkpt == 3:
            logger.info('\nStop at check point 3: BO is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- max_select_bo
        if 0 < self.rin.max_select_bo <= self.n_selection:
            logger.info(f'\nReached max_select_bo: {self.rin.max_select_bo}')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- BO
        backup_cryspy()
        bo_data = (self.init_dscrpt_data, self.opt_dscrpt_data,
                   self.bo_mean, self.bo_var, self.bo_score)
        bo_id_data = (self.n_selection, self.id_queueing,
                      self.id_running, self.id_select_hist)
        from ..BO import bo_next_select
        bo_next_select.next_select(self.rin, self.rslt_data,
                                   bo_id_data, bo_data, noprint)

    def next_select_LAQA(self):
        # ---------- flag for next selection or generation
        if not self.go_next_sg:
            logger.info('\nLAQA is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- check point 3
        if self.rin.stop_chkpt == 3:
            logger.info('\nStop at check point 3: LAQA is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- selection of LAQA
        backup_cryspy()
        laqa_id_data = (self.id_queueing, self.id_running,
                        self.id_select_hist)
        laqa_data = (self.tot_step_select, self.laqa_step, self.laqa_struc,
                     self.laqa_energy, self.laqa_bias, self.laqa_score)
        from ..LAQA import laqa_next_selection
        laqa_next_selection.next_selection(self.rin, laqa_id_data, laqa_data)


#
#  End Ctrl_job class
#


def prepare_jobfile(rin, cid, work_path):
    if not os.path.isfile('./calc_in/' + rin.jobfile):
        logger.error('Could not find ./calc_in' + rin.jobfile)
        raise SystemExit(1)
    with open('./calc_in/' + rin.jobfile, 'r') as f:
        lines = f.readlines()
    lines2 = []
    for line in lines:
        lines2.append(line.replace('CrySPY_ID', str(cid)))
    with open(work_path + rin.jobfile, 'w') as f:
        f.writelines(lines2)


def submit_next_struc(rin, cid, work_path, wait=False):
    # ---------- submit job
    os.chdir(work_path)    # cd work_path
    # ---------- submit
    with open('sublog', 'w') as logf:
        if not wait:
            subprocess.Popen([rin.jobcmd, rin.jobfile],
                             stdout=logf, stderr=logf)
        else:    # wait
            subprocess.run([rin.jobcmd, rin.jobfile],
                           stdout=logf, stderr=logf)
    # ---------- write stat_job
    with open('stat_job', 'w') as f:
        f.write(f'{cid:<6}    # Structure ID\n')
        f.write(f'{1:<6}    # Stage\n')
        f.write('submitted\n')
    os.chdir('../../')    # go back to csp root dir


def regist_opt(
        rin,
        cid,
        work_path,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        opt_struc,
        energy,
        magmom,
        check_opt,
        ef=None,
        nat=None,
        n_selection=None,
        gen=None
    ):
    '''
    Common part in ctrl_collect_*
    '''
    # ---------- get initial spg info
    try:
        spg_sym, spg_num = init_struc_data[
            cid].get_space_group_info(symprec=rin.symprec)
    except TypeError:
        spg_num = 0
        spg_sym = None

    # ---------- success
    if opt_struc is not None:
        # ------ get opt spg info
        try:
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(
                symprec=rin.symprec)
        except:
            spg_num_opt = 0
            spg_sym_opt = None
        # ------ out opt_struc
        out_poscar({cid:opt_struc}, './data/opt_POSCARS')
        try:
            out_cif(opt_struc, cid, work_path,
                    './data/opt_CIFS.cif', rin.symprec)
        except TypeError:
            logger.warning('failed to write opt_CIF')
    # ---------- error
    else:
        spg_num_opt = 0
        spg_sym_opt = None

    # ---------- register opt_struc
    opt_struc_data[cid] = opt_struc
    pkl_data.save_opt_struc(opt_struc_data)

    # ---------- register rslt
    if rin.algo in ['RS', 'LAQA']:
        rslt_data.loc[cid] = [spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo == 'BO':
        rslt_data.loc[cid] = [n_selection,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo in ['EA']:
        rslt_data.loc[cid] = [gen,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo in ['EA-vc']:
        rslt_data.loc[cid] = [gen,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, ef, nat, magmom, check_opt]
    pkl_data.save_rslt(rslt_data)
    if rin.algo != 'EA-vc':
        out_rslt(rslt_data)
    else:    # EA-vc
        out_rslt(rslt_data, order_ef=True)

    # ---------- return
    return opt_struc_data, rslt_data


def update_status(cid, id_queueing, id_running, operation):

    # ---------- read stat
    stat = io_stat.stat_read()

    # ---------- update status
    if operation == 'submit':
        id_running.append(cid)
        id_queueing.remove(cid)
        io_stat.set_stage(stat, cid, 1)
    elif operation == 'fin':
        if cid in id_queueing:
            id_queueing.remove(cid)
        if cid in id_running:
            id_running.remove(cid)
        io_stat.clean_id(stat, cid)
    else:
        logger.error('operation is wrong')
        raise SystemExit(1)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save id_queueing and id_running
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)

    # ---------- return
    return id_queueing, id_running


def mv_fin(cid):
    if not os.path.isdir(f'work/fin/{cid}'):
        shutil.move(f'work/{cid}', 'work/fin/')
    else:    # rename for recalc
        for i in itertools.count(1):
            if not os.path.isdir(f'work/fin/{cid}_{i}'):
                shutil.move(f'work/{cid}',
                            f'work/fin/{cid}_{i}')
                break


def next_gen_EA(
        rin,
        gen,
        go_next_sg,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data=None,
        struc_mol_id=None,
    ):

    # ---------- log
    logger.info(f'Done generation {gen}')

    # ---------- EA-vc: calc convex hull
    if rin.algo == 'EA-vc':
        hdist_data = pkl_data.load_hdist_data()
        pd_data = pkl_data.load_pd_data()
        if gen not in hdist_data:
            logger.info(f'Calculate convex hull for generation {gen}')
            from ..EA.calc_hull import calc_convex_hull
            pd, hdist = calc_convex_hull(
                atype=rin.atype,
                gen=gen,
                end_point=rin.end_point,
                rslt_data=rslt_data,
                nat_data=nat_data,
                show_max=rin.show_max,
                label_stable=rin.label_stable,
                vmax=rin.vmax,
                bottom_margin=rin.bottom_margin,
                emax_ea=rin.emax_ea,
                emin_ea=rin.emin_ea,
            )
            out_hdist(gen, hdist, nat_data)
            hdist_data[gen] = hdist
            pd_data[gen] = pd
            pkl_data.save_hdist_data(hdist_data)
            pkl_data.save_pd_data(pd_data)

    # ---------- flag for next selection or generation
    if not go_next_sg:
        logger.info('\nEA is ready')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- check point 3
    if rin.stop_chkpt == 3:
        logger.info('\nStop at check point 3: EA is ready')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- maxgen_ea
    if 0 < rin.maxgen_ea <= gen:
        logger.info(f'\nReached maxgen_ea: {rin.maxgen_ea}')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- EA
    backup_cryspy()
    from ..EA import ea_next_gen
    ea_next_gen.next_gen(
        rin,
        gen,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        nat_data,
        struc_mol_id,
    )