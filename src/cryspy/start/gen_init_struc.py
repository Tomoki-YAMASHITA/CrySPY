from logging import getLogger
import os

from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, get_mol_data, out_poscar

# import later
#from ..RS.gen_struc_RS import gen_pyxtal
#from ..RS.gen_struc_RS import random_generation


logger = getLogger('cryspy')


def gen_init_struc(rin, prev_nstruc, comm, mpi_rank, mpi_size):
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
    nstruc_list, offset_list = _divide_task(rin.tot_struc, mpi_size, prev_nstruc)

    # ---------- log
    if mpi_rank == 0:
        logger.info('# ------ generate structures')

    # ---------- EA-vc
    vc = True if rin.algo == 'EA-vc' else False

    # ---------- pyxtal
    if not (rin.spgnum == 0 or rin.use_find_wy):
        from ..RS.gen_struc_RS import gen_pyxtal
        # ------ crystal
        if rin.struc_mode == 'crystal':
            init_struc_data = gen_pyxtal.gen_struc(
                                  rin,
                                  nstruc=nstruc_list[mpi_rank],
                                  mindist=mindist,
                                  id_offset=offset_list[mpi_rank],
                                  vc=vc
                              )
            struc_mol_id = {}     # not used, just for return
        # ------ molecular crystal
        elif rin.struc_mode == 'mol':
            mol_data = get_mol_data(rin.mol_file)
            init_struc_data, struc_mol_id = gen_pyxtal.gen_struc_mol(
                                                rin,
                                                nstruc=nstruc_list[mpi_rank],
                                                mindist=mindist,
                                                mol_data=mol_data,
                                                id_offset=offset_list[mpi_rank]
                                            )
        # ------ molecular crystal breaking symmetry
        elif rin.struc_mode == 'mol_bs':
            mol_data = get_mol_data(rin.mol_file)
            init_struc_data, struc_mol_id = gen_pyxtal.gen_struc_mol_break_sym(
                                                rin,
                                                nstruc=nstruc_list[mpi_rank],
                                                mindist=mindist,
                                                mindist_dummy=mindist_dummy,
                                                mol_data=mol_data,
                                                id_offset=offset_list[mpi_rank]
                                            )

    # ---------- w/o pyxtal
    else:
        from ..RS.gen_struc_RS import random_generation
        if rin.spgnum == 0:
            init_struc_data = random_generation.gen_wo_spg(
                                  rin,
                                  nstruc=nstruc_list[mpi_rank],
                                  mindist=mindist,
                                  id_offset=offset_list[mpi_rank],
                                  vc=vc
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
                                  rin,
                                  nstruc=nstruc_list[mpi_rank],
                                  mindist=mindist,
                                  id_offset=offset_list[mpi_rank],
                                  fwpath=fwpath,
                                  mpi_rank=mpi_rank,
                                  vc=vc
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

    # ------ write init_POSCARS
    if mpi_rank == 0:
        out_poscar(init_struc_data, './data/init_POSCARS')

    # ---------- return
    # only init_struc_data in rank0 is important
    return init_struc_data, struc_mol_id


def _divide_task(n, size, init_n = 0):
    '''
    For example, n = 10, size = 4,
    we would like to divide task as follows:
        rank0 = [0, 1, 2] <-- 3
        rank1 = [3, 4, 5] <-- 3
        rank2 = [6, 7]    <-- 2
        rank3 = [8, 9]    <-- 2

    ntask_list = [3, 3, 2, 2]
    offset_list = [0, 3, 6, 8]

    if size == 1:
        ntask_list = [n]
        offset_list = [0]

    case: n = 20, size = 4, init_n = 10
        rank0 = [10, 11, 12] <-- 3
        rank1 = [13, 14, 15] <-- 3
        rank2 = [16, 17]    <-- 2
        rank3 = [18, 19]    <-- 2

    ntask_list = [3, 3, 2, 2]
    offset_list = [10, 13, 16, 18]

    if size == 1:
        ntask_list = [10] <-- n - init_n
        offset_list = [10] <-- init_n
    '''
    q = (n - init_n)//size
    r = (n - init_n)%size
    ntask_list = []
    offset_list = []
    offset = init_n
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