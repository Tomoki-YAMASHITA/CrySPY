#!/usr/bin/env python3
import argparse
from pathlib import Path
import sqlite3

from cryspy.high_throughput.db.record import (
    Status,
    reset_result,
    select_status,
)
from cryspy.high_throughput.db.sqlite import connect_db


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
        'dbfile',
        nargs='?',
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
        with connect_db(str(db_path)) as conn:
            for record_id in args.ids:
                current_status = select_status(conn, record_id)
                if current_status is None:
                    parser.error(f'ID {record_id} does not exist')
                if current_status == Status.RUNNING:
                    parser.error(f'ID {record_id} is RUNNING')

            for record_id in args.ids:
                current_status = select_status(conn, record_id)
                reset_result(conn, record_id, next_status)
                print(
                    f'ID {record_id}: '
                    f'{current_status.name} -> {next_status.name}'
                )

            conn.commit()

    except sqlite3.Error as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()