'''
Restart CrySPY high-throughput mode
'''

from logging import getLogger

import numpy as np

from ...IO import diff_input, pkl_data
from ...IO.read_input import ReadInput
from ...RS.rs_gen import gen_random
from ..db.record import (
    count_records,
    insert_init_struc,
    reset_running_struc,
)
from ..db.sqlite import connect_db
from .check_input import check_input


logger = getLogger('cryspy')


def restart():
    # ---------- read input
    try:
        rin = ReadInput()
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- check input
    check_input(rin)

    # ---------- check input change
    try:
        pin = pkl_data.load_input()
        diff_input.diff_in(rin, pin)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- count records
    with connect_db() as conn:
        n_records = count_records(conn)

    # ---------- check tot_struc
    if rin.tot_struc < n_records:
        logger.error(
            f'tot_struc < number of records ({n_records})'
        )
        raise SystemExit(1)

    # ---------- append structures
    if n_records < rin.tot_struc:
        logger.info('# ---------- Append structures')

        # ------ RNG
        rng = None
        if rin.seed is not None:
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

            if restored:
                logger.info('Restore RNG state')
            else:
                logger.info('Initialize RNG with seed from input')
                rng = np.random.default_rng(rin.seed)
            logger.info(f'RNG seed: {rin.seed}')

        # ------ generate structures
        nstruc = rin.tot_struc - n_records
        init_struc_data = gen_random(
            rin=rin,
            nstruc=nstruc,
            id_offset=n_records + 1,
            comm=None,
            mpi_rank=0,
            mpi_size=1,
            rng=rng,
        )

        # ------ insert structures
        with connect_db() as conn:
            for struc in init_struc_data.values():
                insert_init_struc(conn, struc, rin.nat)

        # ------ save RNG state
        if rng is not None:
            rng_state_data = (rng.bit_generator.state, rin.seed)
            pkl_data.save_rng_state(rng_state_data)

        logger.info(f'Appended {nstruc} structures')

    # ---------- save input
    pkl_data.save_input(rin)

    # ---------- reset running structures
    with connect_db() as conn:
        nreset = reset_running_struc(conn)
    if nreset > 0:
        logger.warning(
            f'Reset {nreset} running structures to waiting'
        )

    # ---------- return
    return rin