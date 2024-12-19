'''
Restart CrySPY
'''

from datetime import datetime
from logging import getLogger
import os

from ..IO import diff_input, pkl_data
from ..IO.read_input import ReadInput
from ..RS.rs_gen import gen_random
from ..util.utility import get_version, backup_cryspy
from ..util.struc_util import out_poscar

# ---------- import later
#from ..RS import rs_restart
#from ..BO import bo_restart
#from ..LAQA import laqa_restart
#from ..EA import ea_append


logger = getLogger('cryspy')

def restart(comm, mpi_rank, mpi_size):
    if mpi_rank == 0:
        logger.info('\n\n\nRestart CrySPY ' + get_version() + '\n\n')

    # ########## MPI start
    if mpi_size > 1:
        comm.barrier()
    # ---------- read input and check the change
    if mpi_rank == 0:
        logger.info('read input, cryspy.in')
    try:
        # all processes read input data
        rin = ReadInput()    # read input data, cryspy,in
    except Exception as e:
        if mpi_rank == 0:
            logger.error(e)
        raise SystemExit(1)
    if mpi_rank == 0:
        try:
            # only rank0 compares current and previous input
            pin = pkl_data.load_input()    # load previous input data from input_data.pkl
            diff_input.diff_in(rin, pin)           # compare current and previous input
            pkl_data.save_input(rin)       # save input data to input_data.pkl
        except Exception as e:
            if mpi_size > 1:
                comm.Abort(1)
            logger.error(e)
            raise SystemExit(1)

    # ------ load init_struc_data for appending structures
    # In EA, one can not change tot_struc, so struc_mol_id need not be considered here
    # _append_struc is not allowed in EA and EA-vc either
    init_struc_data = pkl_data.load_init_struc()
    append_flag = False

    # ---------- append structures
    if rin.algo not in ['EA', 'EA-vc']:
        if len(init_struc_data) < rin.tot_struc:
            # ------ append_flag
            append_flag = True
            # ------ backup
            if mpi_rank == 0:
                backup_cryspy()
            # ------ barrier for MPI
            if mpi_size > 1:
                comm.barrier()
            # ------ append struc.
            if mpi_rank == 0:
                prev_nstruc = len(init_struc_data)
            #        init_struc_data is saved in _append_struc
            init_struc_data = _append_struc(rin, init_struc_data, comm, mpi_rank, mpi_size)
        elif rin.tot_struc < len(init_struc_data):
            logger.error('tot_struc < len(init_struc_data)')
            raise SystemExit(1)

    # ---------- append structures by EA (option)
    # not support MPI
    if rin.algo not in ['EA', 'EA-vc']:
        if rin.append_struc_ea:
            # ------ append_flag
            append_flag = True
            if mpi_rank == 0:
                # ------ backup
                backup_cryspy()
    #
    # struc_mol_id has not developed yet here
    #if rin.struc_mode in ['mol', 'mol_bs']:
    #    struc_mol_id = pkl_data.load_struc_mol_id()
    #
                from ..EA import ea_append
                prev_nstruc = len(init_struc_data)
                # init_struc_data is saved in ea_append.append_struc()
                init_struc_data = ea_append.append_struc(rin, init_struc_data)

    # ---------- post append
    if append_flag:
        if mpi_rank == 0:
            # ------ RS
            if rin.algo == 'RS':
                from ..RS import rs_restart
                rs_restart.restart(rin, prev_nstruc)
            # ------ BO
            if rin.algo == 'BO':
                from ..BO import bo_restart
                bo_restart.restart(rin, init_struc_data, prev_nstruc)
            # ------ LAQA
            if rin.algo == 'LAQA':
                from ..LAQA import laqa_restart
                laqa_restart.restart(rin, prev_nstruc)
            # ------ remove lock_cryspy
            os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- return
    return rin, init_struc_data


def _append_struc(rin, init_struc_data, comm, mpi_rank, mpi_size):
    # ---------- log
    if mpi_rank == 0:
        logger.info('# ---------- Append structures')
    if mpi_size > 1:
        comm.barrier()

    # ---------- time
    if mpi_rank == 0:
        time_start = datetime.today()

    # ---------- generate structures
    # only init_struc_data in rank0 is important
    nstruc = rin.tot_struc - len(init_struc_data)
    tmp_struc_data, tmp_struc_mol_id = gen_random(rin, nstruc, len(init_struc_data),
                                        comm, mpi_rank, mpi_size)

    if mpi_rank == 0:
        out_poscar(tmp_struc_data, './data/init_POSCARS')
        # ---------- update and save
        init_struc_data.update(tmp_struc_data)
        # ToDo: struc_mol_id
        pkl_data.save_init_struc(init_struc_data)

        # ---------- time
        time_end = datetime.today()
        etime = time_end - time_start
        logger.info(f'Elapsed time for structure generation: {etime}')

    # ---------- return
    # only init_struc_data in rank0 is important
    return init_struc_data