#!/usr/bin/env python3
import argparse
from pathlib import Path
import sqlite3

from cryspy.high_throughput.db.record import (
    Status,
    select_ids_by_status,
)
from cryspy.high_throughput.db.sqlite import connect_db


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'


def format_id_ranges(ids: list[int]) -> str:
    """Format consecutive record IDs as ranges."""

    # ---------- no IDs
    if not ids:
        return 'None'

    # ---------- initialize ranges
    ranges = []
    start = ids[0]
    previous = ids[0]

    # ---------- group consecutive IDs
    for record_id in ids[1:]:
        if record_id == previous + 1:
            previous = record_id
            continue

        if start == previous:
            ranges.append(str(start))
        else:
            ranges.append(f'{start}-{previous}')

        start = record_id
        previous = record_id

    # ---------- final range
    if start == previous:
        ranges.append(str(start))
    else:
        ranges.append(f'{start}-{previous}')

    # ---------- return
    return ', '.join(ranges)


def print_status(
    conn: sqlite3.Connection,
) -> None:
    """Print status."""

    # ---------- select IDs
    ids_by_status = select_ids_by_status(conn)

    # ---------- format
    total = sum(len(ids) for ids in ids_by_status.values())
    labels = ['Total'] + [status.name for status in Status]
    label_width = max(len(label) for label in labels)
    value_width = len(str(total))

    # ---------- print status
    print(f'{"Total":<{label_width}} {total:>{value_width}}')
    for status in Status:
        ids = ids_by_status[status]
        print(
            f'{status.name:<{label_width}} '
            f'{len(ids):>{value_width}}: '
            f'{format_id_ranges(ids)}'
        )


def main():
    """Print the current status of CrySPY high-throughput mode."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Show CrySPY high-throughput status'
    )
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    args = parser.parse_args()

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- print header
    print(f'Database: {db_path}')

    # ---------- print status
    try:
        with connect_db(str(db_path)) as conn:
            print_status(conn)
    except sqlite3.Error as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()