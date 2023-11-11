from logging import getLogger
import os
import sys

from ..IO import read_input as rin
from ..util.utility import check_fwpath
from ..util.struc_util import set_mindist, out_poscar


logger = getLogger('cryspy')

def gen_init_struc(init_struc_data, struc_mol_id, comm, mpi_rank, mpi_size):
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
    mindist = set_mindist(rin.mindist, rin.mindist_factor, False, mpi_rank)

    # ---------- nstruc, offset for MPI
    nstruc_list, offset_list = _divide_task(rin.tot_struc, mpi_size, len(init_struc_data))

    # ---------- log
    if mpi_rank == 0:
        logger.info('# ------ generate structures')

    # ---------- pyxtal
    if not (rin.spgnum == 0 or rin.use_find_wy):
        from ..RS.gen_struc_RS.gen_pyxtal import Rnd_struc_gen_pyxtal
        rsgx = Rnd_struc_gen_pyxtal(mindist=mindist)
        # ------ crystal
        if rin.struc_mode == 'crystal':
            if not rin.algo == 'EA-vc':
                rsgx.gen_struc(nstruc=nstruc_list[mpi_rank], id_offset=offset_list[mpi_rank])
            else:    # vc
                rsgx.gen_struc(nstruc=nstruc_list[mpi_rank], id_offset=offset_list[mpi_rank],
                               vc=True)
        # ------ molecular crystal
        elif rin.struc_mode == 'mol':
            rsgx.set_mol()
            rsgx.gen_struc_mol(nstruc=nstruc_list[mpi_rank], id_offset=offset_list[mpi_rank])
        # ------ molecular crystal breaking symmetry
        elif rin.struc_mode == 'mol_bs':
            if mpi_rank == 0:
                logger.info('# -- mindist_mol_bs')
            mindist_dummy = set_mindist(rin.mindist_mol_bs, rin.mindist_mol_bs_factor,
                                        dummy=True, mpi_rank=mpi_rank)
            rsgx.set_mol()
            rsgx.gen_struc_mol_break_sym(nstruc=nstruc_list[mpi_rank],
                                            mindist_dummy=mindist_dummy,
                                            id_offset=offset_list[mpi_rank])
        # ------ init_struc_data
        if mpi_size > 1:
            # -- parallel
            # only init_struc_data in rank0 is important
            if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
                struc_dict, mol_id_dict = _gather_struc(rsgx.init_struc_data, rsgx.struc_mol_id,
                                                        comm, mpi_rank)
                init_struc_data.update(struc_dict)
                struc_mol_id.update(mol_id_dict)
            else:
                struc_dict, _ = _gather_struc(rsgx.init_struc_data, None, comm, mpi_rank)
                init_struc_data.update(struc_dict)
            # init_POSCARS
            if mpi_rank == 0:
                for cid, struc in struc_dict.items():
                    out_poscar(struc, cid, './data/init_POSCARS')
        else:
            # -- serial
            init_struc_data.update(rsgx.init_struc_data)
            if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
                struc_mol_id.update(rsgx.struc_mol_id)
            # init_POSCARS
            for cid, struc in rsgx.init_struc_data.items():
                out_poscar(struc, cid, './data/init_POSCARS')
    # ---------- w/o pyxtal
    else:
        from ..RS.gen_struc_RS.random_generation import Rnd_struc_gen
        rsg = Rnd_struc_gen(mindist=mindist)
        if rin.spgnum == 0:
            if not rin.algo == 'EA-vc':
                rsg.gen_wo_spg(nstruc=nstruc_list[mpi_rank], id_offset=offset_list[mpi_rank])
            else:    # vc
                rsg.gen_wo_spg(nstruc=nstruc_list[mpi_rank], id_offset=offset_list[mpi_rank], vc=True)
        else:
            # ---- findwy
            # -- check fwpath
            fwpath = check_fwpath(rin.fwpath)
            if mpi_rank == 0:
                os.makedirs('tmp_gen_struc', exist_ok=True)
            if mpi_size > 1:
                comm.barrier()
            if not rin.algo == 'EA-vc':
                rsg.gen_with_find_wy(nstruc=nstruc_list[mpi_rank],
                                     id_offset=offset_list[mpi_rank],
                                     fwpath=fwpath, mpi_rank=mpi_rank)
            else:    # vc
                rsg.gen_with_find_wy(nstruc=nstruc_list[mpi_rank],
                                     id_offset=offset_list[mpi_rank],
                                     fwpath=fwpath, mpi_rank=mpi_rank, vc=True)
        # ------ init_struc_data
        if mpi_size > 1:
            # -- parallel
            struc_dict, _ = _gather_struc(rsg.init_struc_data, None, comm, mpi_rank)
            init_struc_data.update(struc_dict)
            # init_POSCARS
            if mpi_rank == 0:
                for cid, struc in struc_dict.items():
                    out_poscar(struc, cid, './data/init_POSCARS')
        else:
            # -- serial
            init_struc_data.update(rsg.init_struc_data)
            # init_POSCARS
            for cid, struc in rsg.init_struc_data.items():
                out_poscar(struc, cid, './data/init_POSCARS')

    # ---------- return
    # only init_struc_data in rank0 is important
    if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
        return init_struc_data, struc_mol_id
    else:
        return init_struc_data


def _divide_task(n, size, init_n = 0):
    '''
    For example, n = 10, size = 4,
    we would like to divide task as follows:
        rank0 = [0, 1, 2] <-- 3
        rank1 = [3, 4, 5] <-- 3
        rank2 = [6, 7]    <-- 2
        rank3 = [8, 9]    <-- 2

    nstruc_list = [3, 3, 2, 2]
    offset_list = [0, 3, 6, 8]

    if size == 1:
        nstruc_list = [n]
        offset_list = [0]

    case: n = 20, size = 4, init_n = 10
        rank0 = [10, 11, 12] <-- 3
        rank1 = [13, 14, 15] <-- 3
        rank2 = [16, 17]    <-- 2
        rank3 = [18, 19]    <-- 2

    nstruc_list = [3, 3, 2, 2]
    offset_list = [10, 13, 16, 18]

    if size == 1:
        nstruc_list = [10] <-- n - init_n
        offset_list = [10] <-- init_n
    '''
    q = (n - init_n)//size
    r = (n - init_n)%size
    nstruc_list = []
    offset_list = []
    offset = init_n
    for i in range(size):
        offset_list.append(offset)
        if r > 0:
            nstruc_list.append(q + 1)
            offset += q + 1
            r = r -1
        else:
            nstruc_list.append(q)
            offset += q
    return nstruc_list, offset_list


def _gather_struc(struc_dict, mol_id_dict, comm, mpi_rank):
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