#!/usr/bin/env python3
'''
Recalculate convex hulls in high-throughput mode
'''
import argparse
import configparser
from contextlib import closing
from pathlib import Path
import sqlite3

from cryspy.IO import diff_input, pkl_data
from cryspy.IO.read_input import ReadInput
from cryspy.high_throughput.db.ea_db import (
    save_ea_hdist,
    save_ea_phase_diagram,
    select_ea_algorithm,
    select_ea_info,
)
from cryspy.high_throughput.db.record_db import (
    Status,
    select_ids_in_generation_range,
    select_status,
)
from cryspy.high_throughput.db.sqlite_db import connect_db
from cryspy.high_throughput.worker.controller_ea import calculate_hull_ea
from cryspy.util.utility import set_logger


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'


def main():
    """Recalculate and save an EA-vc convex hull."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Recalculate a CrySPY high-throughput convex hull'
    )
    parser.add_argument('gen', type=int, help='generation to recalculate')
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    args = parser.parse_args()
    if args.gen < 1:
        parser.error('gen must be a positive integer')

    # ---------- check database
    db_path = Path(args.dbfile)
    if not db_path.is_file():
        parser.error(f'{db_path} does not exist')

    # ---------- lock
    if Path('lock_cryspy').is_file():
        parser.error('lock_cryspy file exists')

    try:
        # ---------- logger
        set_logger(
            noprint=False,
            logfile='log_cryspy-ht',
            errfile='err_cryspy-ht',
        )

        # ---------- read input
        rin = ReadInput(ht=True)
        diff_input.diff_in(rin, pkl_data.load_input())
        if rin.algo != 'EA-vc':
            raise ValueError('This command supports only EA-vc')

        # ---------- recalculate convex hull
        with closing(connect_db(str(db_path))) as conn:
            conn.execute('BEGIN IMMEDIATE')
            if select_ea_algorithm(conn) != 'EA-vc':
                raise ValueError('The database is not for EA-vc')
            if select_ea_info(conn, args.gen) is None:
                raise ValueError(f'Generation {args.gen} was not found')

            # ------ check unfinished calculations
            record_ids = select_ids_in_generation_range(conn, 1, args.gen)
            unfinished_ids = [
                cid for cid in record_ids
                if select_status(conn, cid) in (
                    Status.GENERATING,
                    Status.WAITING,
                    Status.RUNNING,
                )
            ]
            if unfinished_ids:
                raise ValueError(
                    f'Unfinished calculations through generation {args.gen}: '
                    f'{unfinished_ids}'
                )

            phase_diagram, hdist = calculate_hull_ea(
                conn,
                rin,
                args.gen,
                mpl_draw=True,
            )
            save_ea_hdist(conn, args.gen, hdist)
            save_ea_phase_diagram(conn, args.gen, phase_diagram)
            conn.commit()

        print(f'Database: {db_path}')
        print(f'Recalculated convex hull for generation {args.gen}')

    except (OSError, ValueError, sqlite3.Error, configparser.Error) as e:
        parser.error(str(e))


if __name__ == '__main__':
    # ---------- main
    main()
