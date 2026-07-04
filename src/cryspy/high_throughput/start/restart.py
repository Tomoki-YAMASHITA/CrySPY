'''
Restart CrySPY high-throughput mode
'''

from logging import getLogger

from ...IO import diff_input, pkl_data
from ...IO.read_input import ReadInput
from ..db.record import reset_running_struc, select_record_ids
from ..db.sqlite import connect_db
from ..worker.controller_rs import generate_structures_rs
from .check_input import check_input


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

    # ---------- select record IDs
    with connect_db() as conn:
        existing_ids = set(select_record_ids(conn))

    # ---------- check record IDs
    unexpected_ids = sorted(
        cid
        for cid in existing_ids
        if cid < 1 or rin.tot_struc < cid
    )
    if unexpected_ids:
        logger.error(
            f'Database includes IDs outside '
            f'1 to tot_struc: {unexpected_ids}'
        )
        raise SystemExit(1)

    # ---------- missing IDs
    missing_ids = [
        cid
        for cid in range(1, rin.tot_struc + 1)
        if cid not in existing_ids
    ]
    recovery_ids = [
        cid
        for cid in missing_ids
        if cid <= pin.tot_struc
    ]
    append_ids = [
        cid
        for cid in missing_ids
        if pin.tot_struc < cid
    ]

    # ---------- generate structures
    if missing_ids:
        if recovery_ids and append_ids:
            logger.info(
                '# ---------- Recover and append structures'
            )
        elif recovery_ids:
            logger.info(
                '# ---------- Recover missing structures'
            )
        else:
            logger.info('# ---------- Append structures')

        # ------ RNG
        if rin.seed is not None:
            logger.info('Use ID-based RNG')
            logger.info(f'RNG seed: {rin.seed}')

        # ------ generate and insert structures
        generate_structures_rs(
            rin,
            missing_ids,
        )

        # ------ generation results
        if recovery_ids:
            logger.info(
                f'Recovered {len(recovery_ids)} '
                f'missing structures'
            )
        if append_ids:
            logger.info(
                f'Appended {len(append_ids)} structures'
            )

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