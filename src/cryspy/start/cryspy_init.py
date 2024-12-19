'''
Initialize CrySPY
'''

from datetime import datetime
from logging import getLogger
import os

import pandas as pd
import pkg_resources

from ..IO import pkl_data, io_stat, write_input
from ..IO.read_input import ReadInput
from ..RS.rs_gen import gen_random
from ..util.utility import get_version
from ..util.struc_util import out_poscar

# ---------- import later
#from ..RS import rs_init
#from ..BO import bo_init
#from ..LAQA import laqa_init
#from ..EA import ea_init


logger = getLogger('cryspy')


def initialize(comm, mpi_rank, mpi_size):
    # ---------- start
    if mpi_rank == 0:
        logger.info('\n\n\nStart CrySPY ' + get_version() + '\n\n')
    # ---------- check versions
    if mpi_rank == 0:
        logger.info('# ---------- Library version info')
        logger.info(f'pandas version: {pd.__version__}')
        v_pymat = pkg_resources.get_distribution('pymatgen').version
        v_pyxtal = pkg_resources.get_distribution('pyxtal').version
        logger.info(f'pymatgen version: {v_pymat}')
        logger.info(f'pyxtal version: {v_pyxtal}')
    # ---------- read input
    if mpi_rank == 0:
        logger.info('# ---------- Read input file, cryspy.in')
    # ########## MPI start
    if mpi_size > 1:
        comm.barrier()
    try:
        rin = ReadInput()    # read input data, cryspy,in
    except Exception as e:
        if mpi_rank == 0:
            logger.error(e)
        raise SystemExit(1)
    # ########## MPI end
    if mpi_rank == 0:
        # ---------- make data directory
        os.makedirs('data/pkl_data', exist_ok=True)
        # ---------- write and save input file
        write_input.out_input(rin)    # log
        pkl_data.save_input(rin)      # input_data.pkl

    # ---------- generate initial structures
    if not rin.load_struc_flag:
        if mpi_size > 1:
            comm.barrier()
        if mpi_rank == 0:
            # ------ time
            time_start = datetime.today()
        # ########## MPI start
        # ------ structure generation
        # only init_struc_data in rank0 is important
        if rin.algo in ['EA', 'EA-vc']:
            nstruc = rin.n_pop
        else:
            nstruc = rin.tot_struc
        init_struc_data, struc_mol_id = gen_random(rin, nstruc, 0,
                                                       comm, mpi_rank, mpi_size)
        # ########## MPI end
        if mpi_rank == 0:
            # ------ init_POSCARS
            out_poscar(init_struc_data, './data/init_POSCARS', mode='w')
            # ------ save
            pkl_data.save_init_struc(init_struc_data)
            if rin.algo in ['EA', 'EA-vc'] and rin.struc_mode in ['mol', 'mol_bs']:
                pkl_data.save_struc_mol_id(struc_mol_id)
            # ------ time
            time_end = datetime.today()
            etime = time_end - time_start
            logger.info(f'Elapsed time for structure generation: {etime}')
    else:
        # ------ load initial structure
        if mpi_rank == 0:
            logger.info('# --------- Load initial structure data')
            logger.info('Load ./data/pkl_data/init_struc_data.pkl')
            init_struc_data = pkl_data.load_init_struc()
            # -- check
            if rin.algo not in ['EA', 'EA-vc']:
                if not rin.tot_struc == len(init_struc_data):
                    logger.error(f'rin.tot_struc = {rin.tot_struc},'
                                    f' len(init_struc_data) = {len(init_struc_data)}')
                    raise SystemExit(1)
            # -- init_POSCARS
            out_poscar(init_struc_data, './data/init_POSCARS', mode='w')

    if mpi_rank == 0:
        # ---------- initialize stat
        io_stat.stat_init()

        # ---------- initialize opt_struc_data
        opt_struc_data = {}
        pkl_data.save_opt_struc(opt_struc_data)

        # ---------- initialize rslt_data
        rslt_data = pd.DataFrame(columns=['Spg_num', 'Spg_sym',
                                        'Spg_num_opt', 'Spg_sym_opt',
                                        'E_eV_atom', 'Magmom', 'Opt'])
        rslt_data[['Spg_num', 'Spg_num_opt']] = rslt_data[
                                    ['Spg_num', 'Spg_num_opt']].astype(int)
        pkl_data.save_rslt(rslt_data)

        # ---------- initialize for each algorithm
        if rin.algo == 'RS':
            from ..RS import rs_init
            rs_init.initialize(rin)
        elif rin.algo == 'BO':
            from ..BO import bo_init
            bo_init.initialize(rin, init_struc_data, rslt_data)
        elif rin.algo == 'LAQA':
            from ..LAQA import laqa_init
            laqa_init.initialize(rin)
        elif rin.algo in ['EA', 'EA-vc']:
            from ..EA import ea_init
            ea_init.initialize(rin, init_struc_data, rslt_data)

        # ---------- initialize etc
        if rin.kpt_flag:
            kpt_data = {}
            pkl_data.save_kpt(kpt_data)
        if rin.energy_step_flag:
            energy_step_data = {}
            pkl_data.save_energy_step(energy_step_data)
        if rin.struc_step_flag:
            struc_step_data = {}
            pkl_data.save_struc_step(struc_step_data)
        if rin.force_step_flag:
            force_step_data = {}
            pkl_data.save_force_step(force_step_data)
        if rin.stress_step_flag:
            stress_step_data = {}
            pkl_data.save_stress_step(stress_step_data)
