from contextlib import redirect_stdout
from logging import getLogger
import os
from typing import Callable

import numpy as np
from tqdm.notebook import tqdm

from .opt_ase import opt_ase
from ..start import cryspy_restart
from ..IO import io_stat, pkl_data
from ..job.ctrl_job import regist_opt, mv_fin, next_gen_EA

# ---------- import later
#from ..EA.calc_ef import calc_ef


logger = getLogger('cryspy')


def restart_interact(
        njob: int,
        calculator: Callable,
        optimizer: str,
        symmetry: bool = True,
        fmax: float = 0.01,
        steps: int = 2000,
    ) -> None:
    """
    Restart the interactive CrySPY process.

    Args:
        njob (int): Number of jobs to run.
        calculator (Callable): Calculator function to use.
        optimizer (str): Optimizer to use ('BFGS', 'LBFGS', 'FIRE').
        symmetry (bool, optional): Whether to use symmetry. Default is True.
        fmax (float, optional): Maximum force. Default is 0.01.
        steps (int, optional): Number of steps. Default is 2000.

    Raises:
        ValueError: If the calculation code is not 'ASE' or the algorithm is not supported in interactive mode.
    """
    # ---------- ignore MPI
    comm = None
    mpi_rank = 0
    mpi_size = 1

    # ---------- restart
    rin, init_struc_data = cryspy_restart.restart(comm, mpi_rank, mpi_size)
    if rin.calc_code != 'ASE':
        logger.error('Use ASE for calc_code in interactive mode')
        raise ValueError('Use ASE for calc_code in interactive mode')
    if rin.algo not in ['RS', 'EA', 'EA-vc']:
        logger.error(f'algo = {rin.algo} is not supported in interactive mode')
        raise ValueError(f'algo = {rin.algo} is not supported in interactive mode')
    if rin.nstage > 1:
        logger.warning('nstage is ignored in interactive mode')
    if njob == 0:
        njob = rin.njob

    # ---------- mkdir work/fin
    os.makedirs('work/fin', exist_ok=True)

    # ----------
    # ---------- Ctrl_job
    # ----------

    # ---------- load data
    id_queueing = pkl_data.load_id_queueing()
    id_running = pkl_data.load_id_running()
    init_struc_data = pkl_data.load_init_struc()
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()
    if rin.algo in ['EA', 'EA-vc']:
        gen = pkl_data.load_gen()
    else:
        gen = None
    if rin.algo == 'EA-vc':
        nat_data = pkl_data.load_nat_data()

    # ---------- flag for next selection or generation
    if rin.algo in ['BO', 'LAQA', 'EA', 'EA-vc']:
        if (id_queueing or id_running):
            go_next_sg = False
        else:
            go_next_sg = True

    # ----------  check_job
    # id_running is supposed to be empty list at the beginning in interactive mode
    if not rin.stop_next_struc:
        while len(id_running) < njob and id_queueing:
            id_running.append(id_queueing.pop(0))

    # ---------- structure optimization
    # in interactive mode, tmp_id_running is not used
    # because id_running does not change during the for loop below
    for cid in tqdm(id_running):
        # ------ work path
        os.makedirs(f'work/{cid}', exist_ok=True)
        work_path = f'work/{cid}/'
        # ------ struc data
        struc = init_struc_data[cid]
        # ------ optimize structure
        with open(work_path + 'log.out', 'w') as f:    # work/xx/log.out
            with redirect_stdout(f):
                opt_struc, energy, converged = opt_ase(
                                        work_path,
                                        struc,
                                        calculator,
                                        optimizer,
                                        symmetry,
                                        fmax,
                                        steps,
                                    )    # eV/cell

        # ---------- eV/cell --> eV/atom
        if not np.isnan(energy):
            if rin.algo == 'EA-vc':
                nat = nat_data[cid]
            else:
                nat = rin.nat
            natot = sum(nat)
            energy = energy/float(natot)    # eV/cell --> eV/atom

        # ---------- check_opt
        if converged:
            check_opt = 'done'
        else:
            check_opt = 'not_yet'

        # ---------- check
        if np.isnan(energy):
            opt_struc = None
            check_opt = 'no_file'
        if opt_struc is None:
            energy = np.nan
            check_opt = 'no_file'

        # ----------  magmom and check_opt
        magmom = np.nan    # not implemented

        # ---------- calculate Ef for EA-vc
        if rin.algo == 'EA-vc':
            from ..EA.calc_ef import calc_ef
            ef = calc_ef(energy, nat, rin.end_point)
            regist_nat = nat
        else:
            ef = None
            regist_nat = None

        # ---------- register opt data
        opt_struc_data, rslt_data = regist_opt(
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
            ef,
            nat=regist_nat,
            n_selection=None,
            gen=gen,
        )

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
        if rin.algo in ['EA', 'EA-vc']:
            if rin.algo == 'EA':
                nat_data = None
            next_gen_EA(
                rin,
                gen,
                go_next_sg,
                init_struc_data,
                opt_struc_data,
                rslt_data,
                nat_data,
                struc_mol_id=None,
            )
