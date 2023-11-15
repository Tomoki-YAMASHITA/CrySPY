'''
Restart CrySPY
'''

from datetime import datetime
from logging import getLogger
import os

from .gen_init_struc import gen_init_struc
from ..IO import io_stat, pkl_data
from ..IO import read_input as rin
from ..util.utility import get_version, backup_cryspy

# ---------- import later
#from ..RS import rs_restart
#from ..BO import bo_restart
#from ..LAQA import laqa_restart
#from ..EA import ea_append
#from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
#from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen


logger = getLogger('cryspy')

def restart(comm, mpi_rank, mpi_size):
    if mpi_rank == 0:
        logger.info('\n\n\nRestart CrySPY ' + get_version() + '\n\n')
        # ---------- read stat
        stat = io_stat.stat_read()
    else:
        stat = None

    # ########## MPI start
    if mpi_size > 1:
        comm.barrier()
    # ---------- read input and check the change
    try:
        rin.readin()
    except Exception as e:
        if mpi_rank == 0:
            logger.error(str(e.args[0]))
        raise SystemExit(1)
    if mpi_rank == 0:
        try:
            rin.diffinstat(stat)
        except Exception as e:
            if mpi_size > 1:
                comm.Abort(1)
            logger.error(str(e.args[0]))
            raise SystemExit(1)

    # ------ load init_struc_data for appending structures
    # In EA, one can not change tot_struc, so struc_mol_id need not be considered here
    # _append_struc is not allowed in EA and EA-vc, too
    init_struc_data = pkl_data.load_init_struc()

    # ---------- append structures
    if len(init_struc_data) < rin.tot_struc:
        # ------ backup
        if mpi_rank == 0:
            backup_cryspy()
        # ------ barrier for MPI
        if mpi_size > 1:
            comm.barrier()
        # ------ append struc.
        if mpi_rank == 0:
            prev_nstruc = len(init_struc_data)
        init_struc_data = _append_struc(init_struc_data, comm, mpi_rank, mpi_size)
        # ------ post append
        if mpi_rank == 0:
            # -- RS
            if rin.algo == 'RS':
                from ..RS import rs_restart
                rs_restart.restart(stat, prev_nstruc)
            # -- BO
            if rin.algo == 'BO':
                from ..BO import bo_restart
                bo_restart.restart(init_struc_data, prev_nstruc)
            # -- LAQA
            if rin.algo == 'LAQA':
                from ..LAQA import laqa_restart
                laqa_restart.restart(stat, prev_nstruc)
            os.remove('lock_cryspy')
        raise SystemExit()
    elif rin.tot_struc < len(init_struc_data):
        logger.error('tot_struc < len(init_struc_data)')
        raise SystemExit(1)

    # ---------- append structures by EA (option)
    # not support MPI
    if rin.append_struc_ea:
        if mpi_rank == 0:
            # ------ backup
            backup_cryspy()
    #
    # EA-vc is not compatible with append_struc_ea option
    # struc_mol_id has not developed yet here
    #if rin.struc_mode in ['mol', 'mol_bs']:
    #    struc_mol_id = pkl_data.load_struc_mol_id()
    #
            from ..EA import ea_append
            prev_nstruc = len(init_struc_data)
            init_struc_data = ea_append.append_struc(stat, init_struc_data)
            # ------ RS
            if rin.algo == 'RS':
                from ..RS import rs_restart
                rs_restart.restart(stat, prev_nstruc)
            # ------ BO
            if rin.algo == 'BO':
                from ..BO import bo_restart
                bo_restart.restart(init_struc_data, prev_nstruc)
            # ------ LAQA
            if rin.algo == 'LAQA':
                from ..LAQA import laqa_restart
                laqa_restart.restart(stat, prev_nstruc)
            os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- return
    return stat, init_struc_data


def _append_struc(init_struc_data, comm, mpi_rank, mpi_size):
    # ---------- append initial structures
    if mpi_rank == 0:
        logger.info('# ---------- Append structures')
    if mpi_size > 1:
        comm.barrier()
    # ------ time
    if mpi_rank == 0:
        time_start = datetime.today()

    # ---------- gen_init_struc()
    # only init_struc_data in rank0 is important
    init_struc_data = gen_init_struc(init_struc_data, None,
                                        comm, mpi_rank, mpi_size)

    if mpi_rank == 0:
        # ---------- save
        pkl_data.save_init_struc(init_struc_data)
        # ---------- time
        time_end = datetime.today()
        etime = time_end - time_start
        logger.info(f'Elapsed time for structure generation: {etime}')

    # ---------- return
    # only init_struc_data in rank0 is important
    return init_struc_data