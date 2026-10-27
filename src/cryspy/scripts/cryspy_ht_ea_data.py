#!/usr/bin/env python3
import argparse
from contextlib import closing
import csv
from math import isfinite
from pathlib import Path
import sqlite3

from cryspy.high_throughput.db.convert_db import blob_to_array
from cryspy.high_throughput.db.ea_db import (
    select_ea_algorithm,
    select_ea_elite_rows,
    select_ea_generations,
    select_ea_hdist_rows,
    select_ea_info,
    select_ea_origins,
    select_ea_phase_diagram_status,
)
from cryspy.high_throughput.db.sqlite_db import connect_db
from cryspy.util.utility import print_option_box


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'
DEFAULT_MAX_GENERATIONS = 1
DEFAULT_MAX_ROWS = 100
DEFAULT_OUTPUT_NAMES = {
    'info': 'ea_info.csv',
    'origin': 'ea_origin.csv',
    'elite': 'ea_elite.csv',
    'hdist': 'ea_hdist.csv',
    'pd': 'ea_pd.csv',
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
    'll_nat',
    'ul_nat',
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
EA_ELITE_CSV_HEADERS = (
    'Elite_Gen',
    *EA_ELITE_HEADERS,
)
EA_HDIST_HEADERS = (
    'Gen',
    'Struc_ID',
    'Hdist_eV_atom',
)
EA_HDIST_CSV_HEADERS = (
    'Hull_Gen',
    *EA_HDIST_HEADERS,
)
EA_PD_HEADERS = (
    'Gen',
    'PD',
    'Hdist_count',
)


def format_value(value) -> str:
    """Format database value."""

    # ---------- null
    if value is None:
        return 'None'

    # ---------- array blob
    if isinstance(value, bytes):
        return str(tuple(blob_to_array(value).tolist()))

    # ---------- other values
    return str(value)


def print_options(selected_table: str | None) -> None:
    """Print available options in an ASCII box."""

    # ---------- generation selection
    if selected_table == 'hdist':
        options = [
            ('-g GEN', 'select hull generation (default: latest)'),
            ('--dmin D / --dmax D', 'inclusive distance range (eV/atom)'),
            ('-n / --no-sort', 'sort by ID (default: distance)'),
        ]
    else:
        generation_default = '' if selected_table == 'pd' else ' (default: latest)'
        options = [
            (
                '-g GEN / --gmin GEN / --gmax GEN',
                f'select generations{generation_default}',
            ),
        ]

    # ---------- row selection
    if selected_table == 'pd':
        options.extend([
            ('--head N / --tail N', 'first / last N generations'),
            ('--limit N', f'latest N generations ({DEFAULT_MAX_ROWS}; 0: all)'),
            ('-a / --all', 'all matching generations'),
        ])
    else:
        options.extend([
            ('--head N / --tail N', 'first / last N rows'),
            ('--limit N', f'row limit ({DEFAULT_MAX_ROWS}; 0: all)'),
            ('-a / --all', 'all rows in selected generations'),
        ])

    # ---------- table selection
    if selected_table in (None, 'info', 'elite', 'origin'):
        options.extend([
            ('--info / --elite / --origin', 'select a table'),
            ('--hdist / --pd', 'EA-vc only'),
            ('--hdist -g GEN', 'single hull generation; no --gmin/--gmax'),
        ])

    # ---------- output and help
    csv_description = (
        'write CSV (all matching by default)'
        if selected_table == 'pd'
        else 'write CSV'
    )
    options.extend([
        ('-o / --output [PATH]', csv_description),
        ('-h / --help', 'command help'),
    ])

    # ---------- print box
    print_option_box(options)


def select_rows(
    rows: list[tuple],
    limit: int,
    tail: bool = False,
) -> list[tuple]:
    """Select rows for output."""

    # ---------- all rows
    if limit == 0:
        return rows

    # ---------- last rows
    if tail:
        return rows[-limit:]

    # ---------- first rows
    return rows[:limit]


def print_table(
    headers: tuple[str, ...],
    rows: list[tuple],
    limit: int = DEFAULT_MAX_ROWS,
    tail: bool = False,
) -> None:
    """Print rows as a table."""

    # ---------- no rows
    if not rows:
        print('None')
        return

    # ---------- display limit
    display_rows = select_rows(rows, limit, tail)
    if len(display_rows) < len(rows):
        print(
            f'Showing {"last" if tail else "first"} '
            f'{len(display_rows)} of {len(rows)} rows. '
            f'Use --head N, --tail N, or --all to display more.'
        )

    # ---------- format values
    formatted_rows = [
        [format_value(value) for value in row]
        for row in display_rows
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
                row[14],
                row[15],
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
        )
        for row in rows
    ]
    return EA_INFO_HEADERS, prepared_rows


def collect_ea_data(
    conn: sqlite3.Connection,
    generations: list[int],
    selected_table: str | None = None,
) -> tuple[
    tuple[str, ...],
    list[tuple],
    list[tuple],
    list[tuple],
    list[tuple],
]:
    """Collect EA data."""

    # ---------- initialize
    info_rows = []
    origin_rows = []
    elite_rows = []
    hdist_rows = []

    # ---------- hull distance table
    has_hdist = selected_table in (None, 'hdist') and conn.execute(
        """
        SELECT 1
        FROM sqlite_master
        WHERE type = 'table' AND name = 'ea_hdist'
        """
    ).fetchone() is not None

    # ---------- collect generations
    for generation in generations:
        ea_info = None
        if selected_table in (None, 'info') or has_hdist:
            ea_info = select_ea_info(
                conn,
                generation,
            )
        if selected_table in (None, 'info') and ea_info is not None:
            info_rows.append(ea_info)

        if selected_table in (None, 'origin'):
            origins = select_ea_origins(
                conn,
                generation,
            )
            origin_rows.extend(origins)

        if selected_table in (None, 'elite'):
            elites = select_ea_elite_rows(
                conn,
                generation,
            )
            elite_rows.extend(elites)

        if (
            has_hdist
            and ea_info is not None
            and any(value is not None for value in ea_info[9:12])
        ):
            hdist = select_ea_hdist_rows(conn, generation)
            hdist_rows.extend(hdist)

    # ---------- prepare ea_info
    info_headers, info_rows = prepare_ea_info(info_rows)

    # ---------- return
    return (
        info_headers,
        info_rows,
        origin_rows,
        elite_rows,
        hdist_rows,
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

    # ---------- hull distance generation
    if args.hdist:
        return generations[-1:], False

    # ---------- generation range
    if args.gmin is not None or args.gmax is not None:
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
    hdist_rows: list[tuple],
    limit: int = 0,
    tail: bool = False,
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
            EA_ELITE_CSV_HEADERS,
            elite_rows,
        ),
        'hdist': (
            EA_HDIST_CSV_HEADERS,
            hdist_rows,
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
            output_rows = select_rows(rows, limit, tail)
            write_csv(
                output_path,
                headers,
                output_rows,
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
    output_rows = select_rows(rows, limit, tail)
    write_csv(
        output_path,
        headers,
        output_rows,
    )


def print_ea_info(
    info_headers: tuple[str, ...],
    info_rows: list[tuple],
    limit: int = DEFAULT_MAX_ROWS,
    tail: bool = False,
) -> None:
    """Print EA information in compact tables for terminal output."""

    # ---------- population columns
    is_vc = info_headers == EA_VC_INFO_HEADERS
    population_columns = [
        ('Gen', 'Gen'),
        ('Pop', 'Population'),
        ('Cross', 'Crossover'),
        ('Perm', 'Permutation'),
        ('Strain', 'Strain'),
    ]
    if is_vc:
        population_columns.extend([
            ('Add', 'Addition'),
            ('Elim', 'Elimination'),
            ('Subs', 'Substitution'),
        ])
    population_columns.extend([
        ('Rand', 'Random'),
        ('Elite', 'Elite'),
    ])

    # ---------- selection columns
    selection_columns = [
        ('Gen', 'Gen'),
        ('crs_lat', 'crs_lat'),
        ('slct_func', 'slct_func'),
    ]
    groups = [
        ('EA population', population_columns),
        ('EA selection', selection_columns),
    ]

    # ---------- constraint columns
    if is_vc:
        constraint_columns = [
            ('Gen', 'Gen'),
            ('ll_nat', 'll_nat'),
            ('ul_nat', 'ul_nat'),
            ('min_comp', 'min_comp'),
            ('max_comp', 'max_comp'),
        ]
        groups.append(('EA constraints', constraint_columns))

    # ---------- print groups
    for title, columns in groups:
        headers = tuple(label for label, _source in columns)
        indices = [info_headers.index(source) for _label, source in columns]
        rows = [
            tuple(row[index] for index in indices)
            for row in info_rows
        ]
        print('')
        print(title)
        print_table(headers, rows, limit=limit, tail=tail)


def print_ea_data(
    selected_table: str | None,
    info_headers: tuple[str, ...],
    info_rows: list[tuple],
    origin_rows: list[tuple],
    elite_rows: list[tuple],
    hdist_rows: list[tuple],
    elite_generations: list[int],
    hdist_generation: int,
    limit: int = DEFAULT_MAX_ROWS,
    no_sort: bool = False,
    tail: bool = False,
) -> None:
    """Print selected EA data."""

    # ---------- ea_info
    if selected_table in (None, 'info'):
        print_ea_info(
            info_headers,
            info_rows,
            limit=limit,
            tail=tail,
        )

    # ---------- ea_elite
    if selected_table in (None, 'elite'):
        display_elite_rows = select_rows(elite_rows, limit, tail)
        for generation in elite_generations:
            rows = [
                (struc_generation, cid)
                for elite_generation, struc_generation, cid in display_elite_rows
                if elite_generation == generation
            ]
            if not rows and any(row[0] == generation for row in elite_rows):
                continue

            print('')
            print(f'EA elites for generation {generation}')
            print('Gen: structure generation.')
            print('')
            print_table(
                EA_ELITE_HEADERS,
                rows,
                limit=0,
            )

        if len(display_elite_rows) < len(elite_rows):
            print(
                f'Showing {"last" if tail else "first"} '
                f'{len(display_elite_rows)} of {len(elite_rows)} rows. '
                f'Use --head N, --tail N, or --all to display more.'
            )

    # ---------- ea_origin
    if selected_table in (None, 'origin'):
        print('')
        print('EA origins')
        if selected_table is None:
            print('Details: cryspy-ht-ea-data --origin')
        else:
            print_table(
                EA_ORIGIN_HEADERS,
                origin_rows,
                limit=limit,
                tail=tail,
            )

    # ---------- ea_hdist
    if selected_table == 'hdist' or (
        selected_table is None and info_headers == EA_VC_INFO_HEADERS
    ):
        print('')
        if selected_table is None:
            print('EA hull distances')
            print('Details: cryspy-ht-ea-data --hdist')
        else:
            print(f'EA hull distances at generation {hdist_generation}')
            print(f'Number of hull distances: {len(hdist_rows)}')
            print(f'Sorted by: {"id" if no_sort else "distance"}')
            print('Gen: structure generation.')
            print('')

            # Hull generation is in the heading; each row shows structure generation.
            # Round distances only for terminal output; CSV keeps full precision.
            display_hdist_rows = [
                (struc_generation, cid, f'{distance:.8f}')
                for _hull_generation, struc_generation, cid, distance in hdist_rows
            ]
            print_table(
                EA_HDIST_HEADERS,
                display_hdist_rows,
                limit=limit,
                tail=tail,
            )

    # ---------- ea_phase_diagram
    if selected_table is None and info_headers == EA_VC_INFO_HEADERS:
        print('')
        print('EA phase diagrams')
        print('Details: cryspy-ht-ea-data --pd')


def show_pd_data(
    conn: sqlite3.Connection,
    db_path: Path,
    args,
    row_limit: int,
) -> None:
    """Print or write phase diagram storage status."""

    # ---------- algorithm
    if select_ea_algorithm(conn) != 'EA-vc':
        raise ValueError('--pd supports only EA-vc databases')

    # ---------- generation range
    min_generation = args.gmin
    max_generation = args.gmax
    if args.generation is not None:
        min_generation = args.generation
        max_generation = args.generation

    # ---------- select storage status
    # Default to the latest rows; --head selects the earliest rows instead.
    pd_rows = select_ea_phase_diagram_status(
        conn,
        min_generation=min_generation,
        max_generation=max_generation,
        limit=row_limit,
        tail=args.head is None,
    )

    # ---------- CSV output
    if args.output is not None:
        output_path = (
            default_output_directory(db_path) / DEFAULT_OUTPUT_NAMES['pd']
            if args.output == ''
            else Path(args.output)
        )
        write_csv(output_path, EA_PD_HEADERS, pd_rows)
        return

    # ---------- options
    print_options('pd')

    # ---------- storage status
    print('')
    print('EA phase diagrams')
    print(f'Number of displayed generations: {len(pd_rows)}')
    print('PD: Saved or Missing. Hdist_count: number of saved hull distances.')
    print('')
    # Rows are already limited in SQL; display every selected row.
    print_table(EA_PD_HEADERS, pd_rows, limit=0)


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

    # ---------- row selection
    row_selection = parser.add_mutually_exclusive_group()
    row_selection.add_argument(
        '--head',
        type=int,
        help='show first N rows per table',
    )
    row_selection.add_argument(
        '--tail',
        type=int,
        help='show last N rows per table',
    )
    row_selection.add_argument(
        '-a',
        '--all',
        dest='all_rows',
        action='store_true',
        help='show all rows in the selected generations',
    )

    # ---------- generation selection
    parser.add_argument(
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
    table_selection.add_argument(
        '--hdist',
        action='store_true',
        help='show only ea_hdist (hull distances in eV/atom)',
    )
    table_selection.add_argument(
        '--pd',
        action='store_true',
        help=(
            'show PD storage status '
            f'(default: latest {DEFAULT_MAX_ROWS} generations)'
        ),
    )

    # ---------- hull distance range
    parser.add_argument(
        '--dmin',
        type=float,
        help='minimum hull distance in eV/atom (inclusive; --hdist only)',
    )
    parser.add_argument(
        '--dmax',
        type=float,
        help='maximum hull distance in eV/atom (inclusive; --hdist only)',
    )

    # ---------- sort order
    parser.add_argument(
        '-n',
        '--no-sort',
        action='store_true',
        help='sort hull distances by ID instead of distance',
    )

    # ---------- display limit
    row_selection.add_argument(
        '--limit',
        type=int,
        default=DEFAULT_MAX_ROWS,
        help=(
            'maximum rows per table in terminal output '
            f'(0: all, default: {DEFAULT_MAX_ROWS})'
        ),
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

    # ---------- check display limit
    if args.limit < 0:
        parser.error('--limit must be a non-negative integer')

    # ---------- check row selection
    if args.head is not None and args.head <= 0:
        parser.error('--head must be a positive integer')
    if args.tail is not None and args.tail <= 0:
        parser.error('--tail must be a positive integer')

    # ---------- check generation
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

    # ---------- check hull distance generation
    if args.hdist and (args.gmin is not None or args.gmax is not None):
        parser.error(
            '--hdist cannot be used with --gmin or --gmax; '
            'use -g GEN'
        )

    # ---------- check hull distance range
    if args.dmin is not None or args.dmax is not None:
        if not args.hdist:
            parser.error('--dmin and --dmax require --hdist')
    if args.dmin is not None and not isfinite(args.dmin):
        parser.error('--dmin must be finite')
    if args.dmax is not None and not isfinite(args.dmax):
        parser.error('--dmax must be finite')
    if (
        args.dmin is not None
        and args.dmax is not None
        and args.dmin > args.dmax
    ):
        parser.error('--dmin must be less than or equal to --dmax')

    # ---------- row limit
    if args.head is not None:
        row_limit = args.head
    elif args.tail is not None:
        row_limit = args.tail
    elif args.all_rows or args.output is not None:
        row_limit = 0
    else:
        row_limit = args.limit

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
    elif args.hdist:
        selected_table = 'hdist'
    elif args.pd:
        selected_table = 'pd'
    else:
        selected_table = None

    # ---------- print header
    if args.output is None:
        print(f'Database: {db_path}')

    # ---------- EA data
    try:
        with closing(connect_db(str(db_path))) as conn, conn:
            # ------ phase diagram storage status
            if selected_table == 'pd':
                show_pd_data(conn, db_path, args, row_limit)
                return

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
                hdist_rows,
            ) = collect_ea_data(
                conn,
                generations,
                selected_table=selected_table,
            )

        # ---------- hull distance range
        if args.dmin is not None or args.dmax is not None:
            # Filter before sorting and row selection for both terminal and CSV.
            hdist_rows = [
                row for row in hdist_rows
                if (args.dmin is None or row[3] >= args.dmin)
                and (args.dmax is None or row[3] <= args.dmax)
            ]

        # ---------- hull distance order
        if args.output is not None or selected_table == 'hdist':
            if args.no_sort:
                output_hdist_rows = sorted(
                    hdist_rows,
                    key=lambda row: (row[2], row[0]),
                )
            else:
                output_hdist_rows = sorted(
                    hdist_rows,
                    key=lambda row: (row[3], row[0], row[2]),
                )
        else:
            output_hdist_rows = hdist_rows

        # ---------- options
        if args.output is None:
            print_options(selected_table)

        # ---------- default generation limit
        if truncated:
            print('')
            if DEFAULT_MAX_GENERATIONS == 1:
                print('Showing the latest generation.')
            else:
                print(
                    f'Showing the latest {DEFAULT_MAX_GENERATIONS} generations.'
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
                output_hdist_rows,
                limit=row_limit,
                tail=args.tail is not None,
            )
            return

        # ---------- terminal output
        print_ea_data(
            selected_table,
            info_headers,
            info_rows,
            origin_rows,
            elite_rows,
            output_hdist_rows,
            elite_generations=generations,
            hdist_generation=generations[-1],
            limit=row_limit,
            no_sort=args.no_sort,
            tail=args.tail is not None,
        )

    except (OSError, ValueError, sqlite3.Error) as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()
