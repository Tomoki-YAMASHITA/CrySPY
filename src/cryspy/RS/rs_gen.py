from logging import getLogger
import os

from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, get_mol_data

# import later
# from .gen_struc_RS import gen_pyxtal
# from ..gen_struc_RS import random_generation


logger = getLogger('cryspy')


def gen_random(rin, nstruc, id_offset, comm, mpi_rank, mpi_size):
    # ---------- log: num of MPI processes
    if mpi_rank == 0:
        logger.info('# ---------- Initial structure generation')
    if mpi_size > 1 and mpi_rank == 0:
        logger.info(f'Number of MPI processes: {mpi_size}')
    # ---------- mindist
    if mpi_rank == 0:
        logger.info('# ------ mindist')
    if mpi_size > 1:
        comm.barrier()
    mindist = set_mindist(rin.atype, rin.mindist, rin.mindist_factor, rin.struc_mode,
                          False, None, mpi_rank)
    if rin.struc_mode == 'mol_bs':
        mindist_dummy = set_mindist(
            rin.atype,
            rin.mindist_mol_bs,
            rin.mindist_mol_bs_factor,
            rin.struc_mode,
            dummy=True,
            mol_file=rin.mol_file,
            mpi_rank=mpi_rank
        )

    # ---------- nstruc, offset for MPI
    nstruc_list, offset_list = _divide_task(nstruc, mpi_size, id_offset)

    # ---------- log
    if mpi_rank == 0:
        logger.info('# ------ generate structures')

    # ---------- EA-vc
    vc = True if rin.algo == 'EA-vc' else False

    # ---------- pyxtal
    if not (rin.spgnum == 0 or rin.use_find_wy):
        from .gen_struc_RS import gen_pyxtal
        # ------ crystal
        if rin.struc_mode == 'crystal':
            init_struc_data = gen_pyxtal.gen_struc(
                nstruc=nstruc_list[mpi_rank],
                atype=rin.atype,
                nat=rin.nat,
                mindist=mindist,
                spgnum=rin.spgnum,
                symprec=rin.symprec,
                id_offset=offset_list[mpi_rank],
                vol_factor=rin.vol_factor,
                vol_mu=rin.vol_mu,
                vol_sigma=rin.vol_sigma,
                vc=vc,
                ll_nat=rin.ll_nat,
                ul_nat=rin.ul_nat,
                charge=rin.charge,
            )
            struc_mol_id = {}     # not used, just for return
        # ------ molecular crystal
        elif rin.struc_mode == 'mol':
            mol_data = get_mol_data(rin.mol_file)
            init_struc_data, struc_mol_id = gen_pyxtal.gen_struc_mol(
                nstruc=nstruc_list[mpi_rank],
                atype=rin.atype,
                nat=rin.nat,
                mindist=mindist,
                spgnum=rin.spgnum,
                mol_data=mol_data,
                nmol=rin.nmol,
                symprec=rin.symprec,
                id_offset=offset_list[mpi_rank],
                vol_factor=rin.vol_factor,
                vol_mu=rin.vol_mu,
                vol_sigma=rin.vol_sigma,
                timeout_mol=rin.timeout_mol,
            )
        # ------ molecular crystal breaking symmetry
        elif rin.struc_mode == 'mol_bs':
            mol_data = get_mol_data(rin.mol_file)
            init_struc_data, struc_mol_id = gen_pyxtal.gen_struc_mol_break_sym(
                nstruc=nstruc_list[mpi_rank],
                atype=rin.atype,
                nat=rin.nat,
                mindist=mindist,
                mindist_dummy=mindist_dummy,
                spgnum=rin.spgnum,
                mol_data=mol_data,
                nmol=rin.nmol,
                symprec=rin.symprec,
                id_offset=offset_list[mpi_rank],
                vol_factor=rin.vol_factor,
                vol_mu=rin.vol_mu,
                vol_sigma=rin.vol_sigma,
                rot_mol=rin.rot_mol,
                nrot=rin.nrot,
            )
    # ---------- w/o pyxtal
    else:
        from .gen_struc_RS import random_generation
        if rin.spgnum == 0:
            init_struc_data = random_generation.gen_wo_spg(
                nstruc=nstruc_list[mpi_rank],
                atype=rin.atype,
                nat=rin.nat,
                mindist=mindist,
                spgnum=rin.spgnum,
                minlen=rin.minlen,
                maxlen=rin.maxlen,
                dangle=rin.dangle,
                symprec=rin.symprec,
                maxcnt=rin.maxcnt,
                id_offset=offset_list[mpi_rank],
                vol_mu=rin.vol_mu,
                vol_sigma=rin.vol_sigma,
                vc=vc,
                ll_nat=rin.ll_nat,
                ul_nat=rin.ul_nat,
                charge=rin.charge,
            )
            struc_mol_id = {}     # not used, just for return
        else:
            # ---- findwy
            # -- check fwpath
            fwpath = check_fwpath(rin.fwpath)
            if mpi_rank == 0:
                os.makedirs('tmp_gen_struc', exist_ok=True)
            if mpi_size > 1:
                comm.barrier()
            init_struc_data = random_generation.gen_with_find_wy(
                nstruc=nstruc_list[mpi_rank],
                atype=rin.atype,
                nat=rin.nat,
                mindist=mindist,
                spgnum=rin.spgnum,
                minlen=rin.minlen,
                maxlen=rin.maxlen,
                dangle=rin.dangle,
                symprec=rin.symprec,
                maxcnt=rin.maxcnt,
                id_offset=offset_list[mpi_rank],
                vol_mu=rin.vol_mu,
                vol_sigma=rin.vol_sigma,
                fwpath=fwpath,
                mpi_rank=mpi_rank,
                vc=vc,
                ll_nat=rin.ll_nat,
                ul_nat=rin.ul_nat,
                charge=rin.charge,
            )
            struc_mol_id = {}     # not used, just for return

    # ------ gather init_struc_data for MPI
    if mpi_size > 1:
        # -- parallel
        # only init_struc_data in rank0 is important
        init_struc_data, struc_mol_id = _gather_struc(
            rin,
            init_struc_data,
            struc_mol_id,
            comm, mpi_rank
        )

    # ---------- return
    # only init_struc_data in rank0 is important
    return init_struc_data, struc_mol_id


def _divide_task(ntask, size, offset=0):
    '''
    # ---------- args
    ntask: int, total number of tasks
    size: int, number of processes
    offset: int, offset of task number

    # ---------- return
    ntask_list: list, number of tasks for each process
    offset_list: list, offset of task number for each process

    # ---------- description
    For example, ntask = 10, size = 4,
    we would like to divide the task as follows:
        rank0 = [0, 1, 2] <-- 3
        rank1 = [3, 4, 5] <-- 3
        rank2 = [6, 7]    <-- 2
        rank3 = [8, 9]    <-- 2

    ntask_list = [3, 3, 2, 2]
    offset_list = [0, 3, 6, 8]

    ntask = 4, size = 1, offset = 5
        rank0 = [5, 6, 7, 8] <-- 4
        ntask_list = [4]
        offset_list = [5]

    case: ntask = 20, size = 4, init_n = 10
        rank0 = [10, 11, 12, 13, 14] <-- 5
        rank1 = [15, 16, 17, 18, 19] <-- 5
        rank2 = [20, 21, 22, 23, 24] <-- 5
        rank3 = [25, 26, 27, 28, 29] <-- 5

    ntask_list = [5, 5, 5, 5]
    offset_list = [10, 15, 20, 25]
    '''
    q = (ntask)//size
    r = (ntask)%size
    ntask_list = []
    offset_list = []
    for i in range(size):
        offset_list.append(offset)
        if r > 0:
            ntask_list.append(q + 1)
            offset += q + 1
            r = r -1
        else:
            ntask_list.append(q)
            offset += q
    return ntask_list, offset_list


def _gather_struc(rin, struc_dict, mol_id_dict, comm, mpi_rank):
    '''
    only if mpi_size > 1
    '''
    # ---------- struc_data
    data = struc_dict
    data = comm.gather(data, root=0)
    if mpi_rank == 0:
        for d in data:
            struc_dict.update(d)
    else:
        assert data is None
    # ---------- for mol_id
    if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
        data = mol_id_dict
        data = comm.gather(data, root=0)
        if mpi_rank == 0:
            for d in data:
                mol_id_dict.update(d)
        else:
            assert data is None
    # ---------- return
    return struc_dict, mol_id_dict