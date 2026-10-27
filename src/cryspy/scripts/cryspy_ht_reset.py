#!/usr/bin/env python3
import argparse
from contextlib import closing
from pathlib import Path
import shlex
import sqlite3

from cryspy.high_throughput.db.ea_db import (
    delete_ea_hull_data,
    select_ea_algorithm,
    select_latest_ea_generation,
)
from cryspy.high_throughput.db.record_db import (
    Status,
    reset_result,
    select_ids_by_generation,
    select_status,
    select_struc,
)
from cryspy.high_throughput.db.sqlite_db import connect_db


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'


def main():
    """Reset CrySPY high-throughput results."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Reset CrySPY high-throughput results'
    )
    parser.add_argument('ids', nargs='+', type=int, help='structure IDs')
    parser.add_argument(
        '--omit',
        action='store_true',
        help='set status to OMIT instead of WAITING',
    )
    parser.add_argument(
        '--dbfile',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    args = parser.parse_args()

    # ---------- check lock
    if Path('lock_cryspy').is_file():
        parser.error('lock_cryspy file exists')

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- reset results
    next_status = Status.OMIT if args.omit else Status.WAITING
    try:
        with closing(connect_db(str(db_path))) as conn, conn:
            algo = select_ea_algorithm(conn)
            for record_id in args.ids:
                current_status = select_status(conn, record_id)
                if current_status is None:
                    parser.error(f'ID {record_id} does not exist')
                if current_status == Status.RUNNING:
                    parser.error(f'ID {record_id} is RUNNING')
                if not args.omit and current_status == Status.GENERATING:
                    parser.error(f'ID {record_id} is GENERATING')

                # ------ check initial structure for reset
                if not args.omit:
                    initial_struc = select_struc(conn, record_id, initial=True)
                    if (
                        initial_struc is None
                        or initial_struc['atomic_numbers'] is None
                        or initial_struc['lattice'] is None
                        or initial_struc['frac_coords'] is None
                    ):
                        parser.error(
                            f'ID {record_id} has no initial structure; cannot reset'
                        )

            # ------ check generation for reset
            if not args.omit and algo in ('EA', 'EA-vc'):
                generation = select_latest_ea_generation(conn)
                current_ids = set(select_ids_by_generation(conn, generation))
                invalid_ids = [
                    record_id for record_id in args.ids
                    if record_id not in current_ids
                ]
                if invalid_ids:
                    parser.error(
                        f'Only generation {generation} can be reset: '
                        f'invalid IDs = {invalid_ids}'
                    )

            for record_id in args.ids:
                current_status = select_status(conn, record_id)
                reset_result(conn, record_id, next_status)
                print(
                    f'ID {record_id}: '
                    f'{current_status.name} -> {next_status.name}'
                )

            # ------ invalidate convex hull
            if algo == 'EA-vc':
                generation = select_latest_ea_generation(conn)
                if generation is not None:
                    delete_ea_hull_data(
                        conn,
                        generation,
                    )

            conn.commit()

        # ---------- convex hull reminder
        if args.omit and algo == 'EA-vc':
            if generation is not None:
                print(
                    'EA-vc: cleared saved phase diagram and hull distances '
                    f'for generation {generation}.'
                )
                print(
                    'They will be recalculated automatically after '
                    'optimizations finish when cryspy-ht is resumed.'
                )
            command = ['cryspy-ht-calc-convex-hull', 'GEN']
            if str(db_path) != DEFAULT_DB_PATH:
                command.append(str(db_path))
            print(
                'Saved hull data for earlier generations is retained. '
                'Recalculate affected earlier generations manually if needed:'
            )
            print(f'  {shlex.join(command)}')

    except sqlite3.Error as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()
