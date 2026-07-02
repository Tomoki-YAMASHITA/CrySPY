import sqlite3


def initialize_db() -> None:
    """Initialize SQLite database."""

    with sqlite3.connect("data/db_data/rslt_data.db") as conn:

        conn.execute('PRAGMA journal_mode = DELETE')
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA busy_timeout = 5000")
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS records (
                id                INTEGER PRIMARY KEY AUTOINCREMENT,
                status            INTEGER NOT NULL,
                atomic_numbers    BLOB NOT NULL,
                nat               BLOB NOT NULL,
                init_lattice      BLOB NOT NULL,
                init_frac_coords  BLOB NOT NULL,
                opt_lattice       BLOB,
                opt_frac_coords   BLOB,
                energy            REAL
            )
            """
        )
        conn.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_records_status_id
            ON records(status, id)
            """
        )
        conn.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_records_energy_id
            ON records(energy, id)
            """
        )
        conn.commit()


def connect_db(
    db_path: str = 'data/db_data/rslt_data.db',
) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute('PRAGMA foreign_keys = ON')
    conn.execute('PRAGMA busy_timeout = 5000')
    return conn