import os
import shutil

from logging import getLogger
import numpy as np

from ..interface import select_code
from ..IO import read_input as rin
from ..IO import io_stat, pkl_data
from ..IO.out_results import out_rslt
from ..util.struc_util import out_poscar, out_cif

if rin.algo == 'BO':
    from ..BO.select_descriptor import select_descriptor
    from ..BO import bo_next_select
if rin.algo in ['EA', 'EA-vc']:
    from ..EA import ea_next_gen
if rin.algo == 'EA-vc':
    from ..EA.calc_hull import calc_ef


logger = getLogger('cryspy')

class Ctrl_ext:

    def __init__(self, stat, init_struc_data):
        self.stat = stat
        self.init_struc_data = init_struc_data
        self.opt_struc_data = pkl_data.load_opt_struc()
        self.rslt_data = pkl_data.load_rslt()
        self.recheck = False
        self.logic_next = False
        # ---------- for each algorithm
        if rin.algo == 'RS':
            self.id_queueing, self.id_running = pkl_data.load_rs_id()
        elif rin.algo == 'BO':
            (self.n_selection, self.id_queueing,
             self.id_running, self.id_select_hist) = pkl_data.load_bo_id()
            (self.init_dscrpt_data, self.opt_dscrpt_data,
             self.bo_mean, self.bo_var,
             self.bo_score) = pkl_data.load_bo_data()
        elif rin.rin.algo in ['EA', 'EA-vc']:
            (self.gen, self.id_queueing,
             self.id_running) = pkl_data.load_ea_id()
            if rin.algo == 'EA-vc':
                self.nat_data, self.ratio_data = pkl_data.load_ea_vc_data()
            # do not have to load ea_data here.
            # ea_data is used only in ea_next_gen.py

    def check_job(self):
        # ---------- option: recalc
        if rin.recalc:
            self.set_recalc()
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
            self.init_struc_data[cid].to(fmt='cif', filename='ext/queue/{}.cif'.format(cid))
        # ---------- queue --> running
        self.id_running = self.id_queueing[:]
        self.id_queueing = []
        io_stat.set_id(self.stat, 'id_queueing', self.id_queueing)
        io_stat.write_stat(self.stat)
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
        if rin.algo == 'RS':
            self.ctrl_collect_rs()
        elif rin.algo == 'BO':
            self.ctrl_collect_bo()
        elif rin.rin.algo in ['EA', 'EA-vc']:
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
        ext_opt_struc_data, ext_energy_data = select_code.collect('dummy', 'dummy', 'dummy')
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
        ext_opt_struc_data, ext_energy_data = select_code.collect('dummy', 'dummy', 'dummy')
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
                tmp_dscrpt = select_descriptor({cid: opt_struc})
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
        ext_opt_struc_data, ext_energy_data = select_code.collect('dummy', 'dummy', 'dummy')
        # ---------- register opt_struc
        ext_magmom, ext_check_opt, ext_spg_sym, ext_spg_num, \
            ext_spg_sym_opt, ext_spg_num_opt = self.regist_opt(ext_opt_struc_data)
        # ---------- save rslt
        for cid, opt_struc in ext_opt_struc_data.items():
            if not rin.algo == 'EA-vc':
                self.rslt_data.loc[cid] = [self.gen,
                                           ext_spg_num[cid], ext_spg_sym[cid],
                                           ext_spg_num_opt[cid], ext_spg_sym_opt[cid],
                                           ext_energy_data[cid], ext_magmom[cid], ext_check_opt[cid]]
            else:    # for EA-vc
                ef = calc_ef(ext_energy_data[cid], self.ratio_data[cid], rin.end_point)
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
                    cid].get_space_group_info(symprec=rin.symprec)
            except TypeError:
                ext_spg_num[cid] = 0
                ext_spg_sym[cid] = None
            # ------ success opt
            if ext_opt_struc_data[cid] is not None:
                # -- get opt spg info
                try:
                    ext_spg_sym_opt[cid], ext_spg_num_opt[cid] = opt_struc.get_space_group_info(
                        symprec=rin.symprec)
                except TypeError:
                    ext_spg_num_opt[cid] = 0
                    ext_spg_sym_opt[cid] = None
                # -- out opt_struc
                out_poscar(opt_struc, cid, './data/opt_POSCARS')
                try:
                    out_cif(opt_struc, cid, './data/',
                            './data/opt_CIFS.cif', rin.symprec)
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

    def next_sg(self):
        '''
        next selection or generation
        '''
        if rin.algo == 'BO':
            self.next_select_BO()
        if rin.rin.algo in ['EA', 'EA-vc']:
            self.next_gen_EA()

    def next_select_BO(self):
        # ---------- log
        logger.info(f'\nDone selection {self.n_selection}')
        # ---------- done all structures
        if len(self.rslt_data) == rin.tot_struc:
            logger.info('\nDone all structures!')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- check point 3
        if rin.stop_chkpt == 3:
            logger.info('\nStop at check point 3: BO is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- max_select_bo
        if 0 < rin.max_select_bo <= self.n_selection:
            logger.info(f'\nReached max_select_bo: {rin.max_select_bo}')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- BO
        bo_data = (self.init_dscrpt_data, self.opt_dscrpt_data,
                   self.bo_mean, self.bo_var, self.bo_score)
        bo_id_data = (self.n_selection, self.id_queueing,
                      self.id_running, self.id_select_hist)
        bo_next_select.next_select(self.stat, self.rslt_data,
                                   bo_id_data, bo_data)

    def next_gen_EA(self):
        # ---------- log
        logger.info(f'\nDone generation {self.gen}')
        # ---------- check point 3
        if rin.stop_chkpt == 3:
            logger.info('\nStop at check point 3: EA is ready')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- maxgen_ea
        if 0 < rin.maxgen_ea <= self.gen:
            logger.info(f'\nReached maxgen_ea: {rin.maxgen_ea}')
            os.remove('lock_cryspy')
            raise SystemExit()
        # ---------- EA
        ea_id_data = (self.gen, self.id_queueing, self.id_running)
        ea_next_gen.next_gen(self.stat, self.init_struc_data,
                             self.opt_struc_data, self.rslt_data, ea_id_data)

    def save_id_data(self):
        # ---------- save id_data
        if rin.algo == 'RS':
            rs_id_data = (self.id_queueing, self.id_running)
            pkl_data.save_rs_id(rs_id_data)
        if rin.algo == 'BO':
            bo_id_data = (self.n_selection, self.id_queueing,
                          self.id_running, self.id_select_hist)
            pkl_data.save_bo_id(bo_id_data)
        if rin.rin.algo in ['EA', 'EA-vc']:
            ea_id_data = (self.gen, self.id_queueing, self.id_running)
            pkl_data.save_ea_id(ea_id_data)

    def save_data(self):
        # ---------- save ??_data
        if rin.algo == 'BO':
            bo_data = (self.init_dscrpt_data, self.opt_dscrpt_data,
                       self.bo_mean, self.bo_var, self.bo_score)
            pkl_data.save_bo_data(bo_data)
        # ea_data is used only in ea_next_gen.py
