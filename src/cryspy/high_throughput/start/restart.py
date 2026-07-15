'''
Restart CrySPY high-throughput mode
'''

from logging import getLogger

from ...IO import diff_input, pkl_data
from ...IO.read_input import ReadInput
from ..db.record import reset_running_struc
from ..db.sqlite import connect_db
from .check_input import check_input
from .restart_ea import restart_ea
from .restart_rs import restart_rs


logger = getLogger('cryspy')


def restart():
    # ---------- read input
    try:
        rin = ReadInput(ht=True)
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

    # ---------- reset running structures
    with connect_db() as conn:
        nreset = reset_running_struc(conn)
    if nreset > 0:
        logger.warning(
            f'Reset {nreset} running structures to waiting'
        )

    # ---------- algorithm-specific restart
    if rin.algo == 'RS':
        restart_rs(rin, pin.tot_struc)
    elif rin.algo == 'EA':
        restart_ea(rin)

    # ---------- save input
    pkl_data.save_input(rin)

    # ---------- return
    return rin
