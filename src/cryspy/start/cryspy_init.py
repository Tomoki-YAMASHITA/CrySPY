'''
Initialize CrySPY
'''

from datetime import datetime
from importlib.metadata import version
from logging import getLogger
import os

import numpy as np
import pandas as pd

from ..IO import pkl_data, io_stat, write_input
from ..IO.read_input import ReadInput
from ..RS.rs_gen import gen_random
from ..util.utility import get_version
from ..util.struc_util import out_poscar
from ..util.struc_validation import validate_loaded_structures

# ---------- import later
#from ..RS import rs_init
#from ..BO import bo_init
#from ..LAQA import laqa_init
#from ..EA import ea_init
#from ..util.charge_neutral import prepare_cn_data
#from ..util.struc_util import get_feasible_composition


logger = getLogger('cryspy')


def initialize(comm=None, mpi_rank=0, mpi_size=1):
    # ---------- start
    if mpi_rank == 0:
        logger.info('\n\n\nStart CrySPY ' + get_version() + '\n\n')
    # ---------- check versions
    if mpi_rank == 0:
        logger.info('# ---------- Library version info')
        logger.info(f'ase version: {version("ase")}')
        logger.info(f'pandas version: {version("pandas")}')
        logger.info(f'pymatgen version: {version("pymatgen")}')
        logger.info(f'pyxtal version: {version("pyxtal")}')
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
            os.remove('lock_cryspy')
        if mpi_size > 1:
            comm.Abort(1)      # stop for MPI
        raise SystemExit(1)
    # ########## MPI end

    # ---------- mkdir, write and save input
    if mpi_rank == 0:
        # ---------- make data directory
        os.makedirs('data/pkl_data', exist_ok=True)
        # ---------- write and save input file
        write_input.out_input(rin)    # log
        pkl_data.save_input(rin)      # input_data.pkl

    # ---------- RNG (seed is for serial debug only)
    rng = None
    if mpi_size == 1 and rin.seed is not None:
        logger.info('# ---------- Initialize RNG with seed from input (serial run)')
        rng = np.random.default_rng(rin.seed)
        logger.info(f'RNG seed: {rin.seed}')
    elif mpi_rank == 0 and rin.seed is not None:
        logger.warning('seed is ignored in MPI mode')

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
            logger.info('# ---------- Composition constraints')
            from ..util.struc_util import get_feasible_composition
            from ..util.visual_util import save_composition_window
            save_composition_window(
                atype=rin.atype,
                gen=1,
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
                f'    ./data/convex_hull/composition_window_1.{rin.fig_format}'
            )

    # ---------- init
    init_struc_data = {}

    # ------ number of required structures
    if rin.algo in ['EA', 'EA-vc']:
        n_target = rin.n_pop
    else:
        n_target = rin.tot_struc

    # ------ load initial structures
    if rin.load_struc_flag:
        if mpi_rank == 0:
            logger.info('# ---------- Load initial structure data')
            logger.info('Load ./data/pkl_data/init_struc_data.pkl')
            try:
                init_struc_data = pkl_data.load_init_struc()
                validate_loaded_structures(rin, init_struc_data)
            except ValueError as e:
                logger.error(e)
                os.remove('lock_cryspy')
                if mpi_size > 1:
                    comm.Abort(1)      # stop for MPI
                raise SystemExit(1)
            n_loaded = len(init_struc_data)
        else:
            n_loaded = None
        # -- share number of loaded structures ########## MPI
        if mpi_size > 1:
            n_loaded = comm.bcast(n_loaded, root=0)
    else:
        n_loaded = 0

    # ------ number of structures to generate
    nstruc = n_target - n_loaded

    # ------ generate structures
    if nstruc > 0:
        if mpi_size > 1:
            comm.barrier()
        if mpi_rank == 0:
            time_start = datetime.today()
            logger.info('# ---------- Initial structure generation')
            if n_loaded > 0:
                logger.info(
                    f'Generate {nstruc} random structures '
                    f'to complete the initial population'
                )
        # ########## MPI start
        # only init_struc_data in rank0 is important
        tmp_struc_data = gen_random(
            rin=rin,
            nstruc=nstruc,
            id_offset=n_loaded,
            comm=comm,
            mpi_rank=mpi_rank,
            mpi_size=mpi_size,
            cn_data=None,    # let gen_random load cn_comb_data.pkl on each MPI rank
            feasible_N=None,
            rng=rng
        )
        # ########## MPI end
        if mpi_rank == 0:
            init_struc_data.update(tmp_struc_data)

    # ------ save initial structures
    if mpi_rank == 0:
        out_poscar(init_struc_data, './data/init_POSCARS', mode='w')
        pkl_data.save_init_struc(init_struc_data)
        # ------ time
        if nstruc > 0:
            time_end = datetime.today()
            etime = time_end - time_start
            logger.info(f'Elapsed time for structure generation: {etime}')

    if mpi_rank == 0:
        # ---------- initialize stat
        io_stat.stat_init()

        # ---------- initialize opt_struc_data
        opt_struc_data = {}
        pkl_data.save_opt_struc(opt_struc_data)

        # ---------- initialize rslt_data
        rslt_data = pd.DataFrame(columns=[
            'Spg_num',
            'Spg_sym',
            'Spg_num_opt',
            'Spg_sym_opt',
            'E_eV_atom',
            'Magmom',
            'Opt',
        ])
        rslt_data[['Spg_num', 'Spg_num_opt']] = rslt_data[
            ['Spg_num', 'Spg_num_opt']].astype(int)
        pkl_data.save_rslt(rslt_data)

        # ---------- initialize for each algorithm
        if rin.algo == 'RS':
            from ..RS import rs_init
            rs_init.initialize(rin)
        elif rin.algo == 'BO':
            from ..BO import bo_init
            bo_init.initialize(rin, init_struc_data, rslt_data, rng=rng)
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

        # ---------- save RNG state
        if mpi_size == 1 and rng is not None:
            rng_state_data = (rng.bit_generator.state, rin.seed)
            pkl_data.save_rng_state(rng_state_data)