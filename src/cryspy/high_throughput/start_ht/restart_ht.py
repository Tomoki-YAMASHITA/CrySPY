'''
Restart CrySPY high-throughput mode
'''

from contextlib import closing
from logging import getLogger

from ...IO import diff_input, pkl_data
from ...IO.load_struc import move_loaded_files
from ...IO.read_input import ReadInput
from ..db.ea_db import (
    delete_ea_hull_data,
    select_latest_ea_generation,
)
from ..db.record_db import (
    has_records,
    reset_running_struc,
)
from ..db.sqlite_db import connect_db, initialize_db
from .check_input_ht import check_input_ht
from .initialize_ht import (
    load_initial_structures_ht,
    save_initial_structures_ht,
)
from .restart_ea_ht import restart_ea_ht
from .restart_rs_ht import restart_rs_ht


logger = getLogger('cryspy')


def restart_ht():
    # ---------- read input
    try:
        rin = ReadInput(ht=True)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- check input
    check_input_ht(rin)

    # ---------- check input change
    try:
        pin = pkl_data.load_input()
        diff_input.diff_in(rin, pin)
    except Exception as e:
        logger.error(e)
        raise SystemExit(1)

    # ---------- initialize database
    initialize_db(rin.algo)

    # ---------- prepare database for restart
    with closing(connect_db()) as conn, conn:
        # ------ reset running structures
        records_exist = has_records(conn)
        nreset = reset_running_struc(conn)

        # ------ invalidate convex hull after energy cutoff changes
        if rin.algo == 'EA-vc' and (
            rin.emin_ea != pin.emin_ea
            or rin.emax_ea != pin.emax_ea
        ):
            generation = select_latest_ea_generation(conn)
            if generation is not None:
                delete_ea_hull_data(conn, generation)
                logger.info(
                    'Energy cutoffs changed: cleared saved phase diagram '
                    f'and hull distances for generation {generation}'
                )
    if nreset > 0:
        logger.warning(
            f'Reset {nreset} running structures to waiting'
        )

    # ---------- recover loaded initial structures
    # An empty database may remain if initialization stopped before saving.
    if not records_exist:
        init_struc_data, loaded_files = load_initial_structures_ht(rin)
        save_initial_structures_ht(rin, init_struc_data)

        # ------ move loaded files
        if loaded_files:
            move_loaded_files(loaded_files)

    # ---------- algorithm-specific restart
    if rin.algo == 'RS':
        restart_rs_ht(rin, pin.tot_struc)
    elif rin.algo in ('EA', 'EA-vc'):
        restart_ea_ht(rin)

    # ---------- save input
    pkl_data.save_input(rin)

    # ---------- return
    return rin
