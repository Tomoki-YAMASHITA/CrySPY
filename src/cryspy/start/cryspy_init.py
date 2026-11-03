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
from ..IO.load_struc import load_structures, move_loaded_files
from ..IO.read_input import ReadInput
from ..RS.rs_gen import gen_random_batch
from ..util.struc_util import out_poscar

# ---------- import later
#from ..RS import rs_init
#from ..BO import bo_init
#from ..LAQA import laqa_init
#from ..EA import ea_init
#from ..util.charge_neutral import prepare_cn_data
#from ..util.struc_util import get_feasible_composition


logger = getLogger('cryspy')


def initialize():
    # ---------- check versions
    logger.info('# ---------- Library version info')
    logger.info(f'ase version: {version("ase")}')
    logger.info(f'pandas version: {version("pandas")}')
    logger.info(f'pymatgen version: {version("pymatgen")}')
    logger.info(f'pyxtal version: {version("pyxtal")}')

    # ---------- read input
    logger.info('# ---------- Read input file, cryspy.in')
    try:
        rin = ReadInput()    # read input data, cryspy,in
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- make data directory
    os.makedirs('data/pkl_data', exist_ok=True)

    # ---------- write and save input file
    write_input.out_input(rin)    # log
    pkl_data.save_input(rin)      # input_data.pkl

    # ---------- RNG
    rng = None
    if rin.seed is not None:
        logger.info('# ---------- Initialize RNG with seed from input')
        rng = np.random.default_rng(rin.seed)
        logger.info(f'RNG seed: {rin.seed}')

    # ---------- vc: prepare charge-neutral data
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
            raise SystemExit(1)

    # ---------- vc: plot composition window
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
    try:
        init_struc_data, loaded_files = load_structures(rin)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)
    n_loaded = len(init_struc_data)

    # ------ number of structures to generate
    nstruc = n_target - n_loaded

    # ------ generate structures
    if nstruc > 0:
        time_start = datetime.today()
        logger.info('# ---------- Initial structure generation')
        if n_loaded > 0:
            logger.info(
                f'Generate {nstruc} random structures '
                f'to complete the initial population'
            )
        tmp_struc_data = gen_random_batch(
            rin=rin,
            nstruc=nstruc,
            id_offset=n_loaded + 1,
            cn_data=None,
            feasible_N=None,
            rng=rng
        )
        init_struc_data.update(tmp_struc_data)

    # ------ save initial structures
    out_poscar(init_struc_data, './data/init_POSCARS', mode='w')
    pkl_data.save_init_struc(init_struc_data)
    # ------ time
    if nstruc > 0:
        time_end = datetime.today()
        etime = time_end - time_start
        logger.info(f'Elapsed time for structure generation: {etime}')

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
    if rng is not None:
        rng_state_data = (rng.bit_generator.state, rin.seed)
        pkl_data.save_rng_state(rng_state_data)

    # ---------- move loaded files
    if loaded_files:
        move_loaded_files(loaded_files)
