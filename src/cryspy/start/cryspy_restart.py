'''
Restart CrySPY
'''

from datetime import datetime
from logging import getLogger

import numpy as np

from ..IO import diff_input, pkl_data
from ..IO.read_input import ReadInput
from ..RS.rs_gen import gen_random_batch
from ..util.utility import backup_cryspy
from ..util.struc_util import out_poscar

# ---------- import later
#from ..RS import rs_restart
#from ..BO import bo_restart
#from ..LAQA import laqa_restart
#from ..EA import ea_append
#from ..util.charge_neutral import prepare_cn_data
#from ..util.struc_util import get_feasible_composition


logger = getLogger('cryspy')


def restart():
    # ---------- read input and check the change
    logger.info('# ---------- Read input file, cryspy.in')
    try:
        rin = ReadInput()    # read input data, cryspy,in
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)
    try:
        pin = pkl_data.load_input()    # load previous input data from input_data.pkl
        diff_input.diff_in(rin, pin)           # compare current and previous input
        pkl_data.save_input(rin)       # save input data to input_data.pkl
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- RNG
    rng = None
    if rin.seed is not None:
        # ------ restore RNG state
        restored = False
        if rin.seed == pin.seed:
            try:
                rng_state, saved_seed = pkl_data.load_rng_state()
            except FileNotFoundError:
                pass
            else:
                if saved_seed == rin.seed:
                    rng = np.random.default_rng(saved_seed)
                    rng.bit_generator.state = rng_state
                    restored = True

        # ------ restore or initialize
        if restored:
            logger.info('# ---------- Restore RNG state')
        else:
            logger.info('# ---------- Initialize RNG with seed from input')
            rng = np.random.default_rng(rin.seed)
            rng_state_data = (rng.bit_generator.state, rin.seed)
            pkl_data.save_rng_state(rng_state_data)

        # ------ log
        logger.info(f'RNG seed: {rin.seed}')

    # ------ load init_struc_data for appending structures
    # _append_struc is not allowed in EA and EA-vc either
    init_struc_data = pkl_data.load_init_struc()
    append_flag = False

    # ---------- append structures
    if rin.algo not in ['EA', 'EA-vc']:
        if len(init_struc_data) < rin.tot_struc:
            if rin.append_struc_ea:
                logger.error(
                    'Increasing tot_struc and append_struc_ea '
                    'cannot be used together'
                )
                raise SystemExit(1)
            # ------ append_flag
            append_flag = True
            # ------ backup
            backup_cryspy()
            # ------ append struc.
            prev_nstruc = len(init_struc_data)
            #        init_struc_data is saved in _append_struc
            init_struc_data = _append_struc(rin, init_struc_data, rng=rng)
        elif rin.tot_struc < len(init_struc_data):
            logger.error('tot_struc < len(init_struc_data)')
            raise SystemExit(1)

    # ---------- append structures by EA (option)
    if rin.algo not in ['EA', 'EA-vc']:
        if rin.append_struc_ea:
            # ------ append_flag
            append_flag = True
            # ------ backup
            backup_cryspy()
            # ------ append struc
            from ..EA import ea_append
            prev_nstruc = len(init_struc_data)
            # init_struc_data is saved in ea_append.append_struc()
            init_struc_data = ea_append.append_struc(rin, init_struc_data, rng)

    # ---------- post append
    if append_flag:
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
        # ------ save RNG state
        if rng is not None:
            rng_state_data = (rng.bit_generator.state, rin.seed)
            pkl_data.save_rng_state(rng_state_data)
        raise SystemExit()

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


def _append_struc(rin, init_struc_data, rng=None):
    # ---------- log
    logger.info('# ---------- Append structures')

    # ---------- time
    time_start = datetime.today()

    # ---------- generate structures
    nstruc = rin.tot_struc - len(init_struc_data)
    tmp_struc_data = gen_random_batch(
        rin=rin,
        nstruc=nstruc,
        id_offset=len(init_struc_data) + 1,
        cn_data=None,
        feasible_N=None,
        rng=rng,
    )

    out_poscar(tmp_struc_data, './data/init_POSCARS')
    # ---------- update and save
    init_struc_data.update(tmp_struc_data)
    pkl_data.save_init_struc(init_struc_data)

    # ---------- time
    time_end = datetime.today()
    etime = time_end - time_start
    logger.info(f'Elapsed time for structure generation: {etime}')

    # ---------- return
    return init_struc_data
