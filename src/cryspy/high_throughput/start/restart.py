'''
Restart CrySPY high-throughput mode
'''

from logging import getLogger

from ...IO.read_input import ReadInput
from ...util.utility import get_version
from ..db.record import reset_running_struc
from ..db.sqlite import connect_db
from .check_input import check_input


logger = getLogger('cryspy')


def restart():
    # ---------- start
    logger.info(
        '\n\n\nRestart CrySPY high-throughput mode '
        + get_version()
        + '\n\n'
    )

    # ---------- read input
    rin = ReadInput()

    # ---------- check input
    check_input(rin)

    # ---------- reset running structures
    with connect_db() as conn:
        nreset = reset_running_struc(conn)
    if nreset > 0:
        logger.warning(
            f'Reset {nreset} running structures to waiting'
        )

    # ---------- return
    return rin