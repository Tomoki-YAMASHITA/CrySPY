'''
Restart CrySPY
'''

from datetime import datetime
from logging import getLogger
import os

import numpy as np

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
#from ..util.charge_neutral import prepare_cn_data
#from ..util.struc_util import get_feasible_composition


logger = getLogger('cryspy')


def restart(comm=None, mpi_rank=0, mpi_size=1):
    # ---------- start
    if mpi_rank == 0:
        logger.info('\n\n\nRestart CrySPY ' + get_version() + '\n\n')

    # ---------- read input and check the change
    # ########## MPI start
    if mpi_size > 1:
        comm.barrier()
    if mpi_rank == 0:
        logger.info('# ---------- Read input file, cryspy.in')
    try:
        # all processes read input data
        rin = ReadInput()    # read input data, cryspy,in
    except Exception as e:
        if mpi_rank == 0:
            logger.error(e)
            os.remove('lock_cryspy')
        if mpi_size > 1:
            comm.Abort(1)      # stop for MPI
        raise SystemExit(1)
    if mpi_rank == 0:
        try:
            # only rank0 compares current and previous input
            pin = pkl_data.load_input()    # load previous input data from input_data.pkl
            diff_input.diff_in(rin, pin)           # compare current and previous input
            pkl_data.save_input(rin)       # save input data to input_data.pkl
        except Exception as e:
            logger.error(e)
            os.remove('lock_cryspy')
            if mpi_size > 1:
                comm.Abort(1)      # stop for MPI
            raise SystemExit(1)    # stop for sereial

    # ---------- RNG (seed is for serial debug only)
    rng = None
    if mpi_size == 1 and rin.seed is not None:
        logger.info('# ---------- Initialize RNG with seed from input (serial run)')
        rng = np.random.default_rng(rin.seed)
        logger.info(f'RNG seed: {rin.seed}')

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
            init_struc_data = _append_struc(rin, init_struc_data, comm, mpi_rank, mpi_size, rng)
        elif rin.tot_struc < len(init_struc_data):
            if mpi_rank == 0:
                logger.error('tot_struc < len(init_struc_data)')
                os.remove('lock_cryspy')
            if mpi_size > 1:
                comm.Abort(1)      # stop for MPI
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
                init_struc_data = ea_append.append_struc(rin, init_struc_data, rng)

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

    # ---------- vc: prepare charge-neutral data
    if mpi_rank == 0:
        if rin.algo == 'EA-vc' and rin.charge is not None:
            try:
                from ..util.charge_neutral import prepare_cn_data
                cn_data = prepare_cn_data(
                    rin.ll_nat,
                    rin.ul_nat,
                    rin.charge,
                    cn_mode=rin.cn_mode,
                    max_cn_grid_points=rin.max_cn_grid_points,
                )
                if cn_data['mode'] == 'enumerate' and len(cn_data['cn_comb']) == 0:
                    raise ValueError('No charge neutral combinations found.')
            except Exception as e:
                logger.error(e)
                os.remove('lock_cryspy')
                if mpi_size > 1:
                    comm.Abort(1)      # stop for MPI
                raise SystemExit(1)    # stop for serial

    # ---------- vc: plot composition window
    if mpi_rank == 0:
        if (
                rin.algo == 'EA-vc'
                and (rin.min_comp is not None or rin.max_comp is not None)
                and len(rin.atype) in (2, 3)
            ):
            gen = pkl_data.load_gen()
            comp_changed = (
                rin.min_comp != pin.min_comp
                or rin.max_comp != pin.max_comp
                or rin.axis_order != pin.axis_order
                or rin.fig_format != pin.fig_format
            )
            if comp_changed:
                logger.info('# ---------- Composition constraints')
                from ..util.struc_util import get_feasible_composition
                from ..util.visual_util import save_composition_window
                next_gen = gen + 1
                save_composition_window(
                    atype=rin.atype,
                    gen=next_gen,
                    min_comp=rin.min_comp,
                    max_comp=rin.max_comp,
                    fig_format=rin.fig_format,
                    axis_order=rin.axis_order,
                )
                feasible_comp = get_feasible_composition(rin.min_comp, rin.max_comp)
                logger.info('Feasible composition range:')
                for symbol, (lower, upper) in zip(rin.atype, feasible_comp):
                    logger.info(f'    {symbol}: {lower:.3f} - {upper:.3f}')
                logger.info(
                    f'Composition window saved as '
                    f'    ./data/convex_hull/composition_window_{next_gen}.{rin.fig_format}'
                )

    # ---------- return
    return rin, init_struc_data, rng


def _append_struc(rin, init_struc_data, comm, mpi_rank, mpi_size, rng=None):
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
    tmp_struc_data, tmp_struc_mol_id = gen_random(
        rin=rin,
        nstruc=nstruc,
        id_offset=len(init_struc_data),
        comm=comm,
        mpi_rank=mpi_rank,
        mpi_size=mpi_size,
        cn_data=None,    # let gen_random load cn_comb_data.pkl on each MPI rank
        feasible_N=None,
        rng=rng,
    )

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
