'''
Restart evolutionary algorithm in CrySPY high-throughput mode
'''

from logging import getLogger

from ...IO.read_input import ReadInput
from ..db.ea import (
    delete_ea_elites_after_generation,
    delete_ea_origins_after_generation,
    select_latest_ea_generation,
)
from ..db.record import (
    delete_records_after_generation,
    select_ids_after_generation,
    select_record_ids,
)
from ..db.sqlite import connect_db
from ..worker.controller_ea import initialize_ea


logger = getLogger('cryspy')


def restart_ea(
    rin: ReadInput,
) -> None:
    """Restart evolutionary algorithm."""

    # ---------- EA generation information
    with connect_db() as conn:
        generation = select_latest_ea_generation(conn)

        # ------ initial generation
        if generation is None:
            existing_ids = set(select_record_ids(conn))

        # ------ incomplete next generation
        else:
            incomplete_ids = select_ids_after_generation(
                conn,
                generation,
            )
            if incomplete_ids:
                logger.info(
                    f'Discard incomplete EA generation '
                    f'{generation + 1}: '
                    f'{len(incomplete_ids)} structures'
                )
                delete_ea_origins_after_generation(
                    conn,
                    generation,
                )
                delete_ea_elites_after_generation(
                    conn,
                    generation,
                )
                delete_records_after_generation(
                    conn,
                    generation,
                )
                conn.commit()

    # ---------- recover initial generation
    if generation is None:
        expected_ids = set(range(1, rin.n_pop + 1))

        # ------ check record IDs
        unexpected_ids = sorted(existing_ids - expected_ids)
        if unexpected_ids:
            logger.error(
                f'Database includes IDs outside '
                f'initial EA population: {unexpected_ids}'
            )
            raise SystemExit(1)

        # ------ recover random structures
        logger.info(
            '# ---------- Recover initial EA generation'
        )
        initialize_ea(
            rin,
            sorted(expected_ids),
        )