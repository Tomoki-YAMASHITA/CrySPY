from logging import getLogger
import os

from ..start import cryspy_restart

from ..IO import io_stat, pkl_data
from ..job.ctrl_job import prepare_jobfile, submit_next_struc, regist_opt, mv_fin, next_gen_EA
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

    # ---------- restart
    rin, init_struc_data = cryspy_restart.restart(comm, mpi_rank, mpi_size)
    if rin.nstage > 1:
        logger.warning('nstage is ignored in interactive mode')
    if njob == 0:
        njob = rin.njob

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
    id_queueing = pkl_data.load_id_queueing()
    id_running = pkl_data.load_id_running()
    if rin.algo == 'RS':
        gen = None
    elif rin.algo == 'EA':
        gen = pkl_data.load_gen()

    # ---------- flag for next selection or generation
    if rin.algo in ['BO', 'LAQA', 'EA', 'EA-vc']:
        if (id_queueing or id_running):
            go_next_sg = False
        else:
            go_next_sg = True

    # ----------  check_job
    if not rin.stop_next_struc:
        while len(id_running) < njob and id_queueing:
            id_running.append(id_queueing.pop(0))

    # in interactive mode, tmp_id_running is not used
    # because id_running does not change during the for loop below
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

        # ---------- collect
        opt_struc, energy, magmom, check_opt = \
            select_code.collect(rin, cid, work_path, rin.nat)
        logger.info(f'{cid}    collect results: E = {energy} eV/atom')

        # ---------- register opt data
        opt_struc_data, rslt_data = regist_opt(
            rin, cid, work_path,
            init_struc_data, opt_struc_data, rslt_data,
            opt_struc, energy, magmom, check_opt, ef=None, n_selection=None, gen=gen)

        # ---------- mv work to fin
        mv_fin(cid)

    # ---------- id
    id_running.clear()
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- next selection or generation
    if not (id_queueing or id_running):
        # ---------- RS
        if rin.algo == 'RS':
            logger.info('\nDone all structures!')
        # ---------- EA
        if rin.algo == 'EA':
            next_gen_EA(
                rin,
                id_queueing,
                id_running,
                gen,
                go_next_sg,
                init_struc_data,
                opt_struc_data,
                rslt_data,
                ea_vc_data=None,
                struc_mol_id=None,
            )