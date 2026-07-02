#!/usr/bin/env python3
import argparse
from pathlib import Path
import sqlite3
import sys

from cryspy.high_throughput.db.convert import raw_to_struc
from cryspy.high_throughput.db.record import select_struc
from cryspy.high_throughput.db.sqlite import connect_db


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'


def main():
    """Extract structures from a CrySPY high-throughput database."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description=(
            'Extract structures from a CrySPY '
            'high-throughput database'
        )
    )
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    parser.add_argument(
        '-i',
        '--id',
        dest='ids',
        type=int,
        nargs='+',
        required=True,
        help='structure IDs',
    )
    parser.add_argument(
        '--initial',
        action='store_true',
        help='extract initial structures instead of optimized structures',
    )
    parser.add_argument(
        '-s',
        '--symmetrized',
        action='store_true',
        help='write symmetrized CIF files',
    )
    parser.add_argument(
        '--tolerance',
        type=float,
        default=0.01,
        help='tolerance for symmetrization (default: 0.01)',
    )
    parser.add_argument(
        '-o',
        '--output',
        help='output filename; available only for one ID',
    )
    args = parser.parse_args()

    # ---------- check options
    if any(record_id < 1 for record_id in args.ids):
        parser.error('ID must be a positive integer')
    if args.tolerance <= 0:
        parser.error('--tolerance must be positive')
    if args.output is not None and len(args.ids) != 1:
        parser.error('--output is available only when one ID is specified')

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- extract structures
    try:
        with connect_db(str(db_path)) as conn:
            for record_id in args.ids:
                raw_struc_dict = select_struc(
                    conn,
                    record_id,
                    initial=args.initial,
                )

                # ------ ID not found
                if raw_struc_dict is None:
                    print(
                        f'Warning: ID {record_id} does not exist. Skipped.',
                        file=sys.stderr,
                    )
                    continue

                # ------ optimized structure not registered
                if (
                    raw_struc_dict['lattice'] is None
                    or raw_struc_dict['frac_coords'] is None
                ):
                    print(
                        f'Warning: Optimized structure for '
                        f'ID {record_id} is not registered. Skipped.',
                        file=sys.stderr,
                    )
                    continue

                # ------ convert structure
                structure = raw_to_struc(raw_struc_dict)

                # ------ output filename
                if args.output is not None:
                    output = args.output
                else:
                    prefix = 'init' if args.initial else 'opt'
                    output = f'{prefix}_{record_id}.cif'

                # ------ write CIF
                if args.symmetrized:
                    structure.to(
                        fmt='cif',
                        filename=output,
                        symprec=args.tolerance,
                    )
                    spg_symbol, spg_number = (
                        structure.get_space_group_info(
                            symprec=args.tolerance
                        )
                    )
                    print(f'Output: {output}')
                    print(
                        f'Space group: '
                        f'{spg_symbol} ({spg_number})'
                    )
                else:
                    structure.to(
                        fmt='cif',
                        filename=output,
                    )
                    print(f'Output: {output}')

    except sqlite3.Error as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()