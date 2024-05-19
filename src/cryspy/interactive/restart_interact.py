from logging import getLogger
import os
import subprocess

from ..IO import diff_input, io_stat, pkl_data
from ..IO.read_input import ReadInput
from ..job.ctrl_job import prepare_jobfile, submit_next_struc, regist_opt, mv_fin
from ..util.utility import get_version
#from ..start.gen_init_struc import gen_init_struc

# ---------- import later (after rin = ReadInput())
#from ..interface import select_code


logger = getLogger('cryspy')


def restart_interact(njob: int):
    # ---------- logger etc
    comm = None
    mpi_rank = 0
    mpi_size = 1
    logger.info('\n\n\nRestart CrySPY ' + get_version() + '\n\n')

    # ---------- read input and check the change
    rin = ReadInput()    # read input data, cryspy,in
    if rin.nstage > 1:
        logger.warning('nstage is ignored in interactive mode')
    pin = pkl_data.load_input()    # load previous input data from input_data.pkl
    diff_input.diff_in(rin, pin)           # compare current and previous input
    pkl_data.save_input(rin)       # save input data to input_data.pkl
    if njob == 0:
        njob = rin.njob


    # ---------- load init_struc_data for appending structures
    # In EA, one can not change tot_struc, so struc_mol_id need not be considered here
    # _append_struc is not allowed in EA and EA-vc either
    init_struc_data = pkl_data.load_init_struc()
    append_flag = False

    # ---------- append structures
    if len(init_struc_data) < rin.tot_struc:
        logger.info('append structure: not implemented yet')
        raise SystemExit()
        # ------ append_flag
        #append_flag = True
        # ------ backup
        #backup_cryspy()
        # ------ append
        #prev_nstruc = len(init_struc_data)
        #init_struc_data = _append_struc(rin, init_struc_data, comm, mpi_rank, mpi_size)

    # ---------- post append
    # if append_flag:
    #     if mpi_rank == 0:
    #         # ------ RS
    #         if rin.algo == 'RS':
    #             from ..RS import rs_restart
    #             rs_restart.restart(rin, prev_nstruc)
    #         # ------ BO
    #         if rin.algo == 'BO':
    #             from ..BO import bo_restart
    #             bo_restart.restart(rin, init_struc_data, prev_nstruc)
    #         # ------ LAQA
    #         if rin.algo == 'LAQA':
    #             from ..LAQA import laqa_restart
    #             laqa_restart.restart(rin, prev_nstruc)
    #         # ------ remove lock_cryspy
    #         os.remove('lock_cryspy')
    #     raise SystemExit()

    # ---------- check calc files in ./calc_in
    from ..interface import select_code    # after rin = ReadInput()
    select_code.check_calc_files(rin)

    # ---------- mkdir work/fin
    os.makedirs('work/fin', exist_ok=True)

    # ----------
    # ---------- Ctrl_job
    # ----------

    # ---------- load data
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()
    if rin.algo == 'RS':
        id_queueing, id_running = pkl_data.load_rs_id()

    # ----------  check_job
    if not rin.stop_next_struc:
        while len(id_running) < njob and id_queueing:
            id_running.append(id_queueing.pop(0))

    # in interactive mode, tmp_id_running is not used
    # because id_running does not chage during the for loop below
    for cid in id_running:
        # --------- mkdir
        os.makedirs(f'work/{cid}', exist_ok=True)
        work_path = f'work/{cid}/'

        # ---------- next_struc
        #            interactive mode always starts with next_struc
        select_code.next_struc(rin, init_struc_data[cid], cid, work_path, rin.nat)

        # ---------- prepare jobfile
        prepare_jobfile(rin, cid, work_path)

        # ---------- submit job
        submit_next_struc(rin, cid, work_path, wait=True)

        # def ctrl_collect(self)
        # if rin.algo == 'RS':
        #     self.ctrl_collect_rs()
        # def ctrl_collect_rs()
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(rin, cid, work_path, rin.nat)
        logger.info(f'{cid}    collect results: E = {energy} eV/atom')

        # ---------- register opt data
        opt_struc_data, rslt_data = regist_opt(
            rin, cid, work_path,
            init_struc_data, opt_struc_data, rslt_data,
            opt_struc, energy, magmom, check_opt, ef=None, n_selection=None, gen=None)

        # ---------- mv work to fin
        mv_fin(cid)

    # ---------- id
    id_running.clear()
    rs_id_data = (id_queueing, id_running)
    pkl_data.save_rs_id(rs_id_data)
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)
