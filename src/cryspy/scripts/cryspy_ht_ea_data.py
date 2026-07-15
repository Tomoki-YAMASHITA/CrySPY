#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
import sqlite3

from cryspy.high_throughput.db.convert import blob_to_array
from cryspy.high_throughput.db.ea import (
    select_ea_elite_ids,
    select_ea_generations,
    select_ea_info,
    select_ea_origins,
)
from cryspy.high_throughput.db.sqlite import connect_db


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'
DEFAULT_MAX_GENERATIONS = 100
DEFAULT_OUTPUT_NAMES = {
    'info': 'ea_info.csv',
    'origin': 'ea_origin.csv',
    'elite': 'ea_elite.csv',
}
EA_INFO_HEADERS = (
    'Gen',
    'Population',
    'Crossover',
    'Permutation',
    'Strain',
    'Random',
    'Elite',
    'crs_lat',
    'slct_func',
    'min_comp',
    'max_comp',
)
EA_VC_INFO_HEADERS = (
    'Gen',
    'Population',
    'Crossover',
    'Permutation',
    'Strain',
    'Addition',
    'Elimination',
    'Substitution',
    'Random',
    'Elite',
    'crs_lat',
    'slct_func',
    'min_comp',
    'max_comp',
)
EA_ORIGIN_HEADERS = (
    'Gen',
    'Struc_ID',
    'Operation',
    'Parent1',
    'Parent2',
    'Attempts',
)
EA_ELITE_HEADERS = (
    'Gen',
    'Struc_ID',
)


def format_value(value) -> str:
    """Format database value."""

    # ---------- null
    if value is None:
        return 'None'

    # ---------- array blob
    if isinstance(value, bytes):
        return str(blob_to_array(value).tolist())

    # ---------- other values
    return str(value)


def print_table(
    headers: tuple[str, ...],
    rows: list[tuple],
) -> None:
    """Print rows as a table."""

    # ---------- no rows
    if not rows:
        print('None')
        return

    # ---------- format values
    formatted_rows = [
        [format_value(value) for value in row]
        for row in rows
    ]

    # ---------- column widths
    widths = [
        max(
            len(header),
            max(len(row[index]) for row in formatted_rows),
        )
        for index, header in enumerate(headers)
    ]

    # ---------- header
    print(
        '  '.join(
            f'{header:<{width}}'
            for header, width in zip(headers, widths)
        )
    )

    # ---------- rows
    for row in formatted_rows:
        print(
            '  '.join(
                f'{value:<{width}}'
                for value, width in zip(row, widths)
            )
        )


def write_csv(
    output_path: Path,
    headers: tuple[str, ...],
    rows: list[tuple],
) -> None:
    """Write rows to a CSV file."""

    # ---------- output directory
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------- write CSV
    with output_path.open(
        'w',
        newline='',
        encoding='utf-8',
    ) as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([
                format_value(value)
                for value in row
            ])

    # ---------- output
    print(f'Output: {output_path}')
    print(f'Number of rows: {len(rows)}')


def prepare_ea_info(
    rows: list[tuple],
) -> tuple[tuple[str, ...], list[tuple]]:
    """Prepare EA information for output."""

    # ---------- EA-vc
    is_vc = any(
        any(value is not None for value in row[9:12])
        for row in rows
    )

    # ---------- EA-vc rows
    if is_vc:
        prepared_rows = [
            (
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[9],
                row[10],
                row[11],
                row[5],
                row[6],
                row[7],
                row[8],
                row[12],
                row[13],
            )
            for row in rows
        ]
        return EA_VC_INFO_HEADERS, prepared_rows

    # ---------- fixed-composition EA rows
    prepared_rows = [
        (
            row[0],
            row[1],
            row[2],
            row[3],
            row[4],
            row[5],
            row[6],
            row[7],
            row[8],
            row[12],
            row[13],
        )
        for row in rows
    ]
    return EA_INFO_HEADERS, prepared_rows


def collect_ea_data(
    conn: sqlite3.Connection,
    generations: list[int],
) -> tuple[
    tuple[str, ...],
    list[tuple],
    list[tuple],
    list[tuple],
]:
    """Collect EA data."""

    # ---------- initialize
    info_rows = []
    origin_rows = []
    elite_rows = []

    # ---------- collect generations
    for generation in generations:
        ea_info = select_ea_info(
            conn,
            generation,
        )
        if ea_info is not None:
            info_rows.append(ea_info)

        origins = select_ea_origins(
            conn,
            generation,
        )
        origin_rows.extend(
            (
                row[1],
                row[0],
                row[2],
                row[3],
                row[4],
                row[5],
            )
            for row in origins
        )

        elite_ids = select_ea_elite_ids(
            conn,
            generation,
        )
        elite_rows.extend(
            (
                generation,
                elite_id,
            )
            for elite_id in elite_ids
        )

    # ---------- prepare ea_info
    info_headers, info_rows = prepare_ea_info(info_rows)

    # ---------- return
    return (
        info_headers,
        info_rows,
        origin_rows,
        elite_rows,
    )


def default_output_directory(
    db_path: Path,
) -> Path:
    """Return default CSV output directory."""

    # ---------- default database
    if db_path.resolve() == Path(DEFAULT_DB_PATH).resolve():
        return Path('data')

    # ---------- specified database
    return Path.cwd()


def select_generations(
    conn: sqlite3.Connection,
    args,
) -> tuple[list[int], bool]:
    """Select generations for output."""

    # ---------- specified generation
    if args.generation is not None:
        return [args.generation], False

    # ---------- available generations
    generations = select_ea_generations(
        conn,
        min_generation=args.gmin,
        max_generation=args.gmax,
    )

    # ---------- head
    if args.head is not None:
        return generations[:args.head], False

    # ---------- tail
    if args.tail is not None:
        return generations[-args.tail:], False

    # ---------- all or generation range
    if (
        args.all_generations
        or args.gmin is not None
        or args.gmax is not None
    ):
        return generations, False

    # ---------- CSV output
    if args.output is not None:
        return generations, False

    # ---------- default latest generations
    truncated = len(generations) > DEFAULT_MAX_GENERATIONS
    return generations[-DEFAULT_MAX_GENERATIONS:], truncated


def output_ea_data(
    db_path: Path,
    output: str,
    selected_table: str | None,
    info_headers: tuple[str, ...],
    info_rows: list[tuple],
    origin_rows: list[tuple],
    elite_rows: list[tuple],
) -> None:
    """Write selected EA data to CSV."""

    # ---------- tables
    tables = {
        'info': (
            info_headers,
            info_rows,
        ),
        'origin': (
            EA_ORIGIN_HEADERS,
            origin_rows,
        ),
        'elite': (
            EA_ELITE_HEADERS,
            elite_rows,
        ),
    }

    # ---------- default output directory
    default_directory = default_output_directory(db_path)

    # ---------- all tables
    if selected_table is None:
        output_directory = (
            default_directory
            if output == ''
            else Path(output)
        )
        output_directory.mkdir(
            parents=True,
            exist_ok=True,
        )

        for table_name, (headers, rows) in tables.items():
            output_path = (
                output_directory
                / DEFAULT_OUTPUT_NAMES[table_name]
            )
            write_csv(
                output_path,
                headers,
                rows,
            )
        return

    # ---------- selected table
    output_path = (
        default_directory
        / DEFAULT_OUTPUT_NAMES[selected_table]
        if output == ''
        else Path(output)
    )
    headers, rows = tables[selected_table]
    write_csv(
        output_path,
        headers,
        rows,
    )


def print_ea_data(
    selected_table: str | None,
    info_headers: tuple[str, ...],
    info_rows: list[tuple],
    origin_rows: list[tuple],
    elite_rows: list[tuple],
) -> None:
    """Print selected EA data."""

    # ---------- ea_info
    if selected_table in (None, 'info'):
        print('')
        print('EA info')
        print_table(
            info_headers,
            info_rows,
        )

    # ---------- ea_origin
    if selected_table in (None, 'origin'):
        print('')
        print('EA origins')
        print_table(
            EA_ORIGIN_HEADERS,
            origin_rows,
        )

    # ---------- ea_elite
    if selected_table in (None, 'elite'):
        print('')
        print('EA elites')
        print_table(
            EA_ELITE_HEADERS,
            elite_rows,
        )


def main():
    """Print CrySPY high-throughput EA data."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Show CrySPY high-throughput EA data'
    )
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )

    # ---------- generation selection
    generation_selection = parser.add_mutually_exclusive_group()
    generation_selection.add_argument(
        '--head',
        type=int,
        help='show first N EA generations',
    )
    generation_selection.add_argument(
        '--tail',
        type=int,
        help='show last N EA generations',
    )
    generation_selection.add_argument(
        '-a',
        '--all',
        dest='all_generations',
        action='store_true',
        help='show all EA generations',
    )
    generation_selection.add_argument(
        '-g',
        '--generation',
        type=int,
        help='show specified EA generation',
    )
    parser.add_argument(
        '--gmin',
        type=int,
        help='minimum EA generation',
    )
    parser.add_argument(
        '--gmax',
        type=int,
        help='maximum EA generation',
    )

    # ---------- table selection
    table_selection = parser.add_mutually_exclusive_group()
    table_selection.add_argument(
        '--info',
        action='store_true',
        help='show only ea_info',
    )
    table_selection.add_argument(
        '--origin',
        action='store_true',
        help='show only ea_origin',
    )
    table_selection.add_argument(
        '--elite',
        action='store_true',
        help='show only ea_elite',
    )

    # ---------- CSV output
    parser.add_argument(
        '-o',
        '--output',
        nargs='?',
        const='',
        metavar='PATH',
        help='write EA data to CSV',
    )
    args = parser.parse_args()

    # ---------- check generation
    if args.head is not None and args.head <= 0:
        parser.error('--head must be a positive integer')
    if args.tail is not None and args.tail <= 0:
        parser.error('--tail must be a positive integer')
    if (
        args.generation is not None
        and args.generation <= 0
    ):
        parser.error('--generation must be a positive integer')
    if args.gmin is not None and args.gmin <= 0:
        parser.error('--gmin must be a positive integer')
    if args.gmax is not None and args.gmax <= 0:
        parser.error('--gmax must be a positive integer')
    if (
        args.gmin is not None
        and args.gmax is not None
        and args.gmin > args.gmax
    ):
        parser.error(
            '--gmin must be less than or equal to --gmax'
        )
    if (
        args.generation is not None
        and (
            args.gmin is not None
            or args.gmax is not None
        )
    ):
        parser.error(
            '--generation cannot be used with '
            '--gmin or --gmax'
        )
    if (
        args.all_generations
        and (
            args.gmin is not None
            or args.gmax is not None
        )
    ):
        parser.error(
            '--all cannot be used with '
            '--gmin or --gmax'
        )

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- selected table
    if args.info:
        selected_table = 'info'
    elif args.origin:
        selected_table = 'origin'
    elif args.elite:
        selected_table = 'elite'
    else:
        selected_table = None

    # ---------- print header
    if args.output is None:
        print(f'Database: {db_path}')

    # ---------- EA data
    try:
        with connect_db(str(db_path)) as conn:
            generations, truncated = select_generations(
                conn,
                args,
            )

            # ------ no EA data
            if not generations:
                print('No EA data found')
                return

            # ------ collect EA data
            (
                info_headers,
                info_rows,
                origin_rows,
                elite_rows,
            ) = collect_ea_data(
                conn,
                generations,
            )

        # ---------- default limit
        if truncated:
            print(
                f'Showing latest '
                f'{DEFAULT_MAX_GENERATIONS} generations. '
                f'Use --head N, --tail N, or --all '
                f'to display more.'
            )

        # ---------- CSV output
        if args.output is not None:
            output_ea_data(
                db_path,
                args.output,
                selected_table,
                info_headers,
                info_rows,
                origin_rows,
                elite_rows,
            )
            return

        # ---------- terminal output
        print_ea_data(
            selected_table,
            info_headers,
            info_rows,
            origin_rows,
            elite_rows,
        )

    except (OSError, ValueError, sqlite3.Error) as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()