#!/usr/bin/env python3
import argparse
import csv
from math import isfinite
from pathlib import Path
import sqlite3

from cryspy.high_throughput.db.record import (
    count_results,
    get_min_energy,
    select_results,
)
from cryspy.high_throughput.db.sqlite import connect_db


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'
DEFAULT_MAX_ROWS = 100


def print_results(
    conn: sqlite3.Connection,
    ids: list[int] = None,
    no_sort: bool = False,
    head: int = None,
    tail: int = None,
    all_rows: bool = False,
    emin: float = None,
    emax: float = None,
    ewin: float = None,
    output: str = None,
) -> None:
    """Print or write results."""

    # ---------- options
    sort_by = 'id' if no_sort else 'energy'
    has_energy_filter = (
        emin is not None
        or emax is not None
        or ewin is not None
    )
    if head is not None:
        limit = head
    elif tail is not None:
        limit = tail
    elif (
        ids is not None
        or all_rows
        or has_energy_filter
        or output is not None
    ):
        limit = None
    else:
        limit = DEFAULT_MAX_ROWS

    # ---------- energy window
    if ewin is not None:
        min_energy = get_min_energy(conn)
        if min_energy is None:
            print('No results found')
            return
        emin = min_energy
        emax = min_energy + ewin

    # ---------- select results
    rows = select_results(
        conn,
        ids=ids,
        emin=emin,
        emax=emax,
        sort_by=sort_by,
        limit=limit,
        tail=tail is not None,
    )

    # ---------- no results
    if not rows:
        print('No results found')
        return

    # ---------- write CSV
    if output is not None:
        with open(
            output,
            'w',
            newline='',
            encoding='utf-8',
        ) as f:
            writer = csv.writer(f)
            writer.writerow([
                'id',
                'Spg_num',
                'Spg_sym',
                'Spg_num_opt',
                'Spg_sym_opt',
                'Num_atom',
                'E_eV_atom',
                'status',
            ])
            for (
                record_id,
                status,
                energy,
                init_spg_num,
                init_spg_sym,
                opt_spg_num,
                opt_spg_sym,
                nat,
            ) in rows:
                writer.writerow([
                    record_id,
                    init_spg_num,
                    init_spg_sym,
                    opt_spg_num,
                    opt_spg_sym,
                    nat,
                    '' if energy is None else energy,
                    status.name,
                ])

        print(f'Output: {output}')
        print(f'Number of results: {len(rows)}')
        return

    # ---------- number of results
    n_results = count_results(
        conn,
        ids=ids,
        emin=emin,
        emax=emax,
    )

    # ---------- format
    energy_values = [
        'None' if row[2] is None else f'{row[2]:.8f}'
        for row in rows
    ]
    id_width = max(
        len('ID'),
        max(len(str(row[0])) for row in rows),
    )
    init_spg_num_width = max(
        len('Spg_num'),
        max(len(str(row[3])) for row in rows),
    )
    init_spg_sym_width = max(
        len('Spg_sym'),
        max(len(str(row[4])) for row in rows),
    )
    opt_spg_num_width = max(
        len('Spg_num_opt'),
        max(len(str(row[5])) for row in rows),
    )
    opt_spg_sym_width = max(
        len('Spg_sym_opt'),
        max(len(str(row[6])) for row in rows),
    )
    nat_width = max(
        len('Num_atom'),
        max(len(str(row[7])) for row in rows),
    )
    energy_width = max(
        len('E_eV_atom'),
        max(len(value) for value in energy_values),
    )
    status_width = max(
        len('Status'),
        max(len(row[1].name) for row in rows),
    )

    # ---------- information
    print(f'Number of results: {n_results}')
    print(f'Sorted by: {sort_by}')
    if (
        limit is not None
        and head is None
        and tail is None
        and n_results > limit
    ):
        print(
            f'Showing first {limit} results. '
            f'Use --head N, --tail N, or --all to display more.'
        )
    print('')

    # ---------- header
    print(
        f'{"ID":>{id_width}}  '
        f'{"Spg_num":>{init_spg_num_width}}  '
        f'{"Spg_sym":<{init_spg_sym_width}}  '
        f'{"Spg_num_opt":>{opt_spg_num_width}}  '
        f'{"Spg_sym_opt":<{opt_spg_sym_width}}  '
        f'{"Num_atom":<{nat_width}}  '
        f'{"E_eV_atom":>{energy_width}}  '
        f'{"Status":<{status_width}}'
    )

    # ---------- results
    for row, energy_value in zip(rows, energy_values):
        (
            record_id,
            status,
            _,
            init_spg_num,
            init_spg_sym,
            opt_spg_num,
            opt_spg_sym,
            nat,
        ) = row
        print(
            f'{record_id:>{id_width}}  '
            f'{init_spg_num:>{init_spg_num_width}}  '
            f'{str(init_spg_sym):<{init_spg_sym_width}}  '
            f'{opt_spg_num:>{opt_spg_num_width}}  '
            f'{str(opt_spg_sym):<{opt_spg_sym_width}}  '
            f'{str(nat):<{nat_width}}  '
            f'{energy_value:>{energy_width}}  '
            f'{status.name:<{status_width}}'
        )


def main():
    """Print CrySPY high-throughput results."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Show CrySPY high-throughput results'
    )
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        '--head',
        type=int,
        help='show first N results',
    )
    selection.add_argument(
        '--tail',
        type=int,
        help='show last N results',
    )
    selection.add_argument(
        '-a',
        '--all',
        dest='all_rows',
        action='store_true',
        help='show all results',
    )
    selection.add_argument(
        '-i',
        '--id',
        type=int,
        nargs='+',
        help='show specified structure IDs',
    )
    parser.add_argument(
        '-n',
        '--no-sort',
        action='store_true',
        help='sort by ID instead of energy',
    )
    parser.add_argument(
        '--emin',
        type=float,
        help='minimum energy in eV/atom',
    )
    parser.add_argument(
        '--emax',
        type=float,
        help='maximum energy in eV/atom',
    )
    parser.add_argument(
        '--ewin',
        type=float,
        help='energy window from the minimum energy in eV/atom',
    )
    parser.add_argument(
        '-o',
        '--output',
        help='write results to a CSV file',
    )
    args = parser.parse_args()

    # ---------- check options
    if args.head is not None and args.head <= 0:
        parser.error('--head must be a positive integer')
    if args.tail is not None and args.tail <= 0:
        parser.error('--tail must be a positive integer')
    if args.emin is not None and not isfinite(args.emin):
        parser.error('--emin must be finite')
    if args.emax is not None and not isfinite(args.emax):
        parser.error('--emax must be finite')
    if args.ewin is not None:
        if not isfinite(args.ewin) or args.ewin < 0:
            parser.error('--ewin must be a non-negative finite number')
        if args.emin is not None or args.emax is not None:
            parser.error(
                '--ewin cannot be used with --emin or --emax'
            )
    if (
        args.emin is not None
        and args.emax is not None
        and args.emin > args.emax
    ):
        parser.error('--emin must be less than or equal to --emax')

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- print header
    if args.output is None:
        print(f'Database: {db_path}')

    # ---------- print results
    try:
        with connect_db(str(db_path)) as conn:
            print_results(
                conn,
                ids=args.id,
                no_sort=args.no_sort,
                head=args.head,
                tail=args.tail,
                all_rows=args.all_rows,
                emin=args.emin,
                emax=args.emax,
                ewin=args.ewin,
                output=args.output,
            )
    except (OSError, sqlite3.Error) as e:
        parser.error(str(e))

if __name__ == '__main__':
    main()