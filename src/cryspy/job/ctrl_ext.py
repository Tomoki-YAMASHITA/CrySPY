from logging import getLogger
import os
import shutil

import numpy as np

from ..interface import select_code
from ..IO import io_stat, pkl_data
from ..IO.out_results import out_rslt, out_hdist
from ..util.utility import backup_cryspy
from ..util.struc_util import out_poscar, out_cif

# ---------- import later
#from ..BO.select_descriptor import select_descriptor
#from ..BO import bo_next_select
#from ..EA import ea_next_gen
#from ..EA.calc_ef import calc_ef


logger = getLogger('cryspy')


class Ctrl_ext:

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
        if rin.algo == 'RS':
            pass
        elif rin.algo == 'BO':
            self.n_selection = pkl_data.load_n_selection()
            self.id_select_hist = pkl_data.load_id_select_hist()
            self.init_dscrpt_data = pkl_data.load_init_dscrpt_data()
            self.opt_dscrpt_data = pkl_data.load_opt_dscrpt_data()
            self.bo_mean = pkl_data.load_bo_mean()
            self.bo_var = pkl_data.load_bo_var()
            self.bo_score = pkl_data.load_bo_score()
        elif rin.rin.algo in ['EA', 'EA-vc']:
            self.gen = pkl_data.load_gen()
            if rin.algo == 'EA-vc':
                self.nat_data = pkl_data.load_nat_data()
                self.ratio_data = pkl_data.load_ratio_data()
                self.hdist_data = pkl_data.load_hdist_data()
            # do not have to load ea_data here.
            # ea_data is used only in ea_next_gen.py

    def check_job(self):
        # ---------- option: recalc
        #if self.rin.recalc is not None:
        #    self.set_recalc()
        # ---------- check job status
        try:
            with open('ext/stat_job') as fstat:
                jstat = fstat.readline()    # submitted or done or ...
            if jstat[0:3] == 'sub':
                self.job_stat = 'submitted'
            elif jstat[0:4] == 'done':
                self.job_stat = 'done'
            elif jstat[0:3] == 'out':
                self.job_stat = 'out'
            elif jstat[0:2] == 'no':
                self.job_stat = 'no queue'
            else:
                self.job_stat = 'else'
        except IOError:
            self.job_stat = 'no_file'

    def handle_job(self):
        logger.info('# ---------- job status')
        if self.job_stat == 'submitted':
            logger.info('still queueing or running')
        elif self.job_stat == 'done':
            logger.info('collect data')
            self.ctrl_done()
        elif self.job_stat == 'out':
            logger.info('write queueing structure data in ext/queue/')
            self.out_queue()
        elif self.job_stat == 'no queue':
            pass
        elif self.job_stat == 'else':
            logger.error('Wrong job_stat')
            raise SystemExit(1)
        elif self.job_stat == 'no_file':
            logger.error('Wrong job_stat')
            raise SystemExit(1)
        else:
            logger.error('Unexpected error')
            raise SystemExit(1)

    def out_queue(self):
        # ---------- out cifs
        os.makedirs('ext/queue', exist_ok=False)    # if dir exists, error --> stop
        for cid in self.id_queueing:
            self.init_struc_data[cid].to(fmt='cif', filename=f'ext/queue/{cid}.cif')
        # ---------- queue --> running
        self.id_running = self.id_queueing[:]
        self.id_queueing = []
        stat = io_stat.stat_read()
        io_stat.set_id(stat, 'id_queueing', self.id_queueing)
        io_stat.write_stat(stat)
        self.save_id_data()
        # ---------- stat_job
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('submitted\n')
        # ---------- mkdir calc_data
        os.makedirs('ext/calc_data', exist_ok=True)

    def ctrl_done(self):
        # ---------- collect result
        self.ctrl_collect()

    def ctrl_collect(self):
        # ---------- each algo
        if self.rin.algo == 'RS':
            self.ctrl_collect_rs()
        elif self.rin.algo == 'BO':
            self.ctrl_collect_bo()
        elif self.rin.algo in ['EA', 'EA-vc']:
            self.ctrl_collect_ea()
        else:
            logger.error('Error, algo')
            raise SystemExit(1)
        # ---------- calc_data --> old_calc_data
        if os.path.isdir('ext/old_calc_data'):
            shutil.rmtree('ext/old_calc_data')
        shutil.move('ext/calc_data', 'ext/old_calc_data')
        # ---------- queue --> old_queue
        if os.path.isdir('ext/old_queue'):
            shutil.rmtree('ext/old_queue')
        shutil.move('ext/queue', 'ext/old_queue')
        # ---------- stat_job
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('no queue\n')
        # ---------- update status
        self.id_running = []
        self.save_id_data()

    def ctrl_collect_rs(self):
        # ---------- get opt data
        ext_opt_struc_data, ext_energy_data = select_code.collect(self.rin, 'dummy', 'dummy', 'dummy')
        # ---------- register opt_struc
        ext_magmom, ext_check_opt, ext_spg_sym, ext_spg_num, \
            ext_spg_sym_opt, ext_spg_num_opt = self.regist_opt(ext_opt_struc_data)
        # ---------- save rslt
        for cid, opt_struc in ext_opt_struc_data.items():
            self.rslt_data.loc[cid] = [ext_spg_num[cid], ext_spg_sym[cid],
                                       ext_spg_num_opt[cid], ext_spg_sym_opt[cid],
                                       ext_energy_data[cid], ext_magmom[cid], ext_check_opt[cid]]
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)

    def ctrl_collect_bo(self):
        # ---------- get opt data
        ext_opt_struc_data, ext_energy_data = select_code.collect(self.rin, 'dummy', 'dummy', 'dummy')
        # ---------- register opt_struc
        ext_magmom, ext_check_opt, ext_spg_sym, ext_spg_num, \
            ext_spg_sym_opt, ext_spg_num_opt = self.regist_opt(ext_opt_struc_data)
        # ---------- save rslt
        for cid, opt_struc in ext_opt_struc_data.items():
            self.rslt_data.loc[cid] = [self.n_selection,
                                       ext_spg_num[cid], ext_spg_sym[cid],
                                       ext_spg_num_opt[cid], ext_spg_sym_opt[cid],
                                       ext_energy_data[cid], ext_magmom[cid], ext_check_opt[cid]]
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)
        # ---------- descriptor
        for cid, opt_struc in ext_opt_struc_data.items():
            if opt_struc is not None:
                # ------ calc descriptor for opt structure
                from ..BO.select_descriptor import select_descriptor
                tmp_dscrpt = select_descriptor(self.rin, {cid: opt_struc})
                # ------ update descriptors
                self.opt_dscrpt_data.update(tmp_dscrpt)
                # ---------- error
            else:
                # ------ update descriptors and non_error_id
                self.opt_dscrpt_data[cid] = None
        # ---------- save bo_data
        self.save_data()

    def ctrl_collect_ea(self):
        # ---------- get opt data
        ext_opt_struc_data, ext_energy_data = select_code.collect(self.rin, 'dummy', 'dummy', 'dummy')
        # ---------- register opt_struc
        ext_magmom, ext_check_opt, ext_spg_sym, ext_spg_num, \
            ext_spg_sym_opt, ext_spg_num_opt = self.regist_opt(ext_opt_struc_data)
        # ---------- save rslt
        for cid, opt_struc in ext_opt_struc_data.items():
            if not self.rin.algo == 'EA-vc':
                self.rslt_data.loc[cid] = [self.gen,
                                           ext_spg_num[cid], ext_spg_sym[cid],
                                           ext_spg_num_opt[cid], ext_spg_sym_opt[cid],
                                           ext_energy_data[cid], ext_magmom[cid], ext_check_opt[cid]]
            else:    # for EA-vc
                from ..EA.calc_ef import calc_ef
                ef = calc_ef(ext_energy_data[cid], self.ratio_data[cid],
                             self.rin.end_point)
                self.rslt_data.loc[cid] = [self.gen,
                                           ext_spg_num[cid], ext_spg_sym[cid],
                                           ext_spg_num_opt[cid], ext_spg_sym_opt[cid],
                                           ext_energy_data[cid], ef, ext_magmom[cid], ext_check_opt[cid]]
        pkl_data.save_rslt(self.rslt_data)
        out_rslt(self.rslt_data)

    def regist_opt(self, ext_opt_struc_data):
        '''
        Common part in ctrl_collect_*
        '''
        # ---------- space group etc.
        ext_magmom = {}
        ext_check_opt = {}
        ext_spg_sym = {}
        ext_spg_num = {}
        ext_spg_sym_opt = {}
        ext_spg_num_opt = {}
        for cid, opt_struc in ext_opt_struc_data.items():
            ext_magmom[cid] = np.nan
            ext_check_opt[cid] = 'no_file'
            # ------ get initial spg info
            try:
                ext_spg_sym[cid], ext_spg_num[cid] = self.init_struc_data[
                    cid].get_space_group_info(symprec=self.rin.symprec)
            except TypeError:
                ext_spg_num[cid] = 0
                ext_spg_sym[cid] = None
            # ------ success opt
            if ext_opt_struc_data[cid] is not None:
                # -- get opt spg info
                try:
                    ext_spg_sym_opt[cid], ext_spg_num_opt[cid] = opt_struc.get_space_group_info(
                        symprec=self.rin.symprec)
                except TypeError:
                    ext_spg_num_opt[cid] = 0
                    ext_spg_sym_opt[cid] = None
                # -- out opt_struc
                out_poscar({cid:opt_struc}, './data/opt_POSCARS')
                try:
                    out_cif(opt_struc, cid, './data/',
                            './data/opt_CIFS.cif', self.rin.symprec)
                except TypeError:
                    logger.info('failed to write opt_CIF')
            # ------ error
            else:
                ext_spg_num_opt[cid] = 0
                ext_spg_sym_opt[cid] = None
        # ---------- save opt_struc_data
        self.opt_struc_data.update(ext_opt_struc_data)
        pkl_data.save_opt_struc(self.opt_struc_data)
        # ---------- return
        return ext_magmom, ext_check_opt, ext_spg_sym, ext_spg_num, ext_spg_sym_opt, ext_spg_num_opt

    def next_sg(self, noprint=False):
        '''
        next selection or generation
        '''
        if self.rin.algo == 'BO':
            self.next_select_BO(noprint)
        if self.rin.algo in ['EA', 'EA-vc']:
            self.next_gen_EA()

    def next_select_BO(self, noprint=False):
        # ---------- log
        logger.info(f'Done selection {self.n_selection}')
        # ---------- done all structures
        if len(self.rslt_data) == self.rin.tot_struc:
            logger.info('\nDone all structures!')
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
        bo_data = (self.init_dscrpt_data, self.opt_dscrpt_data,
                   self.bo_mean, self.bo_var, self.bo_score)
        bo_id_data = (self.n_selection, self.id_queueing,
                      self.id_running, self.id_select_hist)
        from ..BO import bo_next_select
        bo_next_select.next_select(self.rin, self.rslt_data,
                                   bo_id_data, bo_data, noprint)

    def next_gen_EA(self):
        # ---------- log
        logger.info(f'Done generation {self.gen}')
        # ---------- check point 3
        if self.rin.stop_chkpt == 3:
            logger.info('\nStop at check point 3: EA is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- maxgen_ea
        if 0 < self.rin.maxgen_ea <= self.gen:
            if self.rin.algo == 'EA-vc':
                # ------ when gen reaches maxgen_ea,  next generation is not created
                #        so convex hull is calculated here
                #        update only convex hull and hdist, not elite_struc and elite_fitness
                c_rslt = self.rslt_data[self.rslt_data['Gen'] == self.gen]
                c_ids = c_rslt.index.values    # current IDs [array]
                ef_all = self.rslt_data['Ef_eV_atom'].to_dict()    # formation energy of all structures
                from ..EA.calc_hull import calc_convex_hull_2d
                hdist = calc_convex_hull_2d(self.rin, self.ratio_data, ef_all, c_ids, self.gen)
                out_hdist(self.gen, hdist, self.ratio_data)
                self.hdist_data[self.gen] = hdist
                pkl_data.save_hdist_data(self.hdist_data)
            logger.info(f'\nReached maxgen_ea: {self.rin.maxgen_ea}')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- EA
        backup_cryspy()
        ea_id_data = (self.gen, self.id_queueing, self.id_running)
        if self.rin.algo == 'EA-vc':
            ea_vc_data = (self.nat_data, self.ratio_data, self.hdist_data)
        else:
            ea_vc_data = None
        from ..EA import ea_next_gen
        ea_next_gen.next_gen(self.rin, self.init_struc_data, None,
                             self.opt_struc_data, self.rslt_data, ea_id_data, ea_vc_data)

    def save_id_data(self):
        # ---------- save id_data
        pkl_data.save_id_queueing(self.id_queueing)
        pkl_data.save_id_running(self.id_running)
        if self.rin.algo == 'RS':
            pass
        if self.rin.algo == 'BO':
            pkl_data.save_n_selection(self.n_selection)
            pkl_data.save_id_select_hist(self.id_select_hist)
        if self.rin.algo in ['EA', 'EA-vc']:
            pkl_data.save_gen(self.gen)


    def save_data(self):
        # ---------- save ??_data
        if self.rin.algo == 'BO':
            pkl_data.save_init_dscrpt_data(self.init_dscrpt_data)
            pkl_data.save_opt_dscrpt_data(self.opt_dscrpt_data)
            pkl_data.save_bo_mean(self.bo_mean)
            pkl_data.save_bo_var(self.bo_var)
            pkl_data.save_bo_score(self.bo_score)
        # ea_data is used only in ea_next_gen.py
