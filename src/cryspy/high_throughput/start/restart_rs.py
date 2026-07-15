'''
Restart random search in CrySPY high-throughput mode
'''

from logging import getLogger

from ...IO.read_input import ReadInput
from ..db.record import (
    Status,
    select_ids_by_status,
    select_record_ids,
)
from ..db.sqlite import connect_db
from ..worker.controller_rs import generate_structures_rs


logger = getLogger('cryspy')


def restart_rs(
    rin: ReadInput,
    previous_tot_struc: int,
) -> None:
    """Restart random search."""

    # ---------- select record IDs
    with connect_db() as conn:
        existing_ids = set(select_record_ids(conn))
        ids_by_status = select_ids_by_status(conn)
        generating_ids = set(
            ids_by_status[Status.GENERATING]
        )

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
    missing_ids = {
        cid
        for cid in range(1, rin.tot_struc + 1)
        if cid not in existing_ids
    }

    # ---------- IDs to generate
    generate_ids = sorted(
        generating_ids | missing_ids
    )
    recovery_ids = [
        cid
        for cid in generate_ids
        if cid <= previous_tot_struc
    ]
    append_ids = [
        cid
        for cid in generate_ids
        if previous_tot_struc < cid
    ]

    # ---------- generate structures
    if generate_ids:
        if recovery_ids and append_ids:
            logger.info(
                '# ---------- Recover and append structures'
            )
        elif recovery_ids:
            logger.info(
                '# ---------- Recover structures'
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
            generate_ids,
        )

        # ------ generation results
        if recovery_ids:
            logger.info(
                f'Recovered {len(recovery_ids)} structures'
            )
        if append_ids:
            logger.info(
                f'Appended {len(append_ids)} structures'
            )