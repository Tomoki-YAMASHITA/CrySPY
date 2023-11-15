import itertools
from logging import getLogger
import os
#import sys
import subprocess
import shutil
#from datetime import datetime

from ..IO import read_input as rin
from ..IO import io_stat, pkl_data
from ..IO.out_results import out_rslt
from ..util.utility import get_version, backup_cryspy
from ..util.struc_util import out_poscar, out_cif
#from ..start.gen_init_struc import gen_init_struc

# ---------- import later (after rin.readin())
#from ..interface import select_code

logger = getLogger('cryspy')


def restart_interact(njob: int):
    # ---------- logger etc
    comm = None
    mpi_rank = 0
    mpi_size = 1
    logger.info('\n\n\nRestart CrySPY ' + get_version() + '\n\n')

    # ---------- read stat
    stat = io_stat.stat_read()

    # ---------- read input and check the change
    rin.readin()
    rin.diffinstat(stat)
    if njob == 0:
        njob = rin.njob
    from ..interface import select_code

    # ---------- load init_struc_data for appending structures
    # In EA, one can not change tot_struc, so struc_mol_id need not be considered here
    init_struc_data = pkl_data.load_init_struc()

    # ---------- append structures
    if len(init_struc_data) < rin.tot_struc:
        logger.info('append structure: not implemented yet')
        raise SystemExit()
        # ------ backup
        #backup_cryspy()
        
        # ------ append
        #prev_nstruc = len(init_struc_data)
        #init_struc_data = _append_struc(init_struc_data, comm, mpi_rank, mpi_size)


    if rin.algo == 'RS':
        id_queueing, id_running = pkl_data.load_rs_id()
    # ----------

    os.makedirs('work/fin', exist_ok=True)

    # ----------   def check_job
    if not rin.stop_next_struc:
        while len(id_running) < njob and id_queueing:
            id_running.append(id_queueing.pop(0))

    for cid in id_running:
        # tmp_id_running = id_running.copy()
        tmp_opt_struc = pkl_data.load_opt_struc()
        rslt_data = pkl_data.load_rslt()

        # ----- mkdir
        if not os.path.isdir('work/{:06}'.format(cid)):
            os.mkdir('work/{:06}'.format(cid))
        work_path = './work/{:06}/'.format(cid)

        select_code.next_struc(init_struc_data[cid], cid, work_path)

        # def prepare_jobfile(self)
        if not os.path.isfile('./calc_in/' + rin.jobfile):
            raise IOError('Could not find ./calc_in' + rin.jobfile)
        with open('./calc_in/' + rin.jobfile, 'r') as f:
            lines = f.readlines()
        lines2 = []
        for line in lines:
            lines2.append(line.replace('CrySPY_ID', str(cid)))
        with open(work_path + rin.jobfile, 'w') as f:
            f.writelines(lines2)
        # -------------------------

        # def submit_next_struc(self)
        os.chdir(work_path)
        with open('sublog', 'w') as logf:
            subprocess.run([rin.jobcmd, rin.jobfile],
                           stdout=logf, stderr=logf)
        with open('stat_job', 'w') as fwstat:
            fwstat.write('{:<6}    # Structure ID\n'.format(cid))
            fwstat.write('{:<6}    # Stage\n'.format(1))
            fwstat.write('done\n')
        os.chdir('../../')    # go back to csp root dir
        # --------------------------

        # def ctrl_collect(self)
        # if rin.algo == 'RS':
        #     self.ctrl_collect_rs()
        # def ctrl_collect_rs()
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(cid, work_path, rin.nat)
        print('{0}    collect results: E = {1} eV/atom'.format(cid, energy))
        # ---------- register opt_struc
        spg_sym, spg_num, spg_sym_opt, spg_num_opt = regist_opt(
            init_struc_data, opt_struc, tmp_opt_struc, cid, work_path)
        # ---------- save rslt
        rslt_data.loc[cid] = [spg_num, spg_sym,
                              spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
        pkl_data.save_rslt(rslt_data)
        out_rslt(rslt_data)
        # -----------------------------

        mv_fin(cid)

    id_running.clear()
    rs_id_data = (id_queueing, id_running)
    pkl_data.save_rs_id(rs_id_data)
    print("finish")


def regist_opt(init_struc_data, opt_struc, tmp_opt_struc, cid, work_path):
    '''
    Common part in ctrl_collect_*
    '''
    # ---------- get initial spg info
    try:
        spg_sym, spg_num = init_struc_data[cid].get_space_group_info(
            symprec=rin.symprec)
    except TypeError:
        spg_num = 0
        spg_sym = None
    # ---------- success
    if opt_struc is not None:
        # ------ get opt spg info
        try:
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(
                symprec=rin.symprec)
        except TypeError:
            spg_num_opt = 0
            spg_sym_opt = None
        # ------ out opt_struc
        out_poscar(opt_struc, cid, './data/opt_POSCARS')
        try:
            out_cif(opt_struc, cid, work_path,
                    './data/opt_CIFS.cif', rin.symprec)
        except TypeError:
            print('failed to write opt_CIF')
    # ---------- error
    else:
        spg_num_opt = 0
        spg_sym_opt = None
    # ---------- register opt_struc
    # self.opt_struc_data[self.current_id] = opt_struc
    tmp_opt_struc[cid] = opt_struc
    pkl_data.save_opt_struc(tmp_opt_struc)
    # ---------- return
    return spg_sym, spg_num, spg_sym_opt, spg_num_opt


def mv_fin(cid):
    # print("mv_fin work!")
    if not os.path.isdir('work/fin/{0:06}'.format(cid)):
        shutil.move('work/{:06}'.format(cid), 'work/fin/')
        # print("true!")
    else:    # rename for recalc
        for i in itertools.count(1):
            if not os.path.isdir('work/fin/{0:06}_{1}'.format(cid, i)):
                shutil.move('work/{:06}'.format(cid),
                            'work/fin/{0:06}_{1}'.format(cid, i))
                break


def update_status(operation, cid, stat, id_queueing, id_running):
    print("update_status work!")
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
        raise ValueError('operation is wrong')
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

