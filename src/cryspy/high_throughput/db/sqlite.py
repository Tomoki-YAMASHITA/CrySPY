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
                generation        INTEGER,
                status            INTEGER NOT NULL,
                atomic_numbers    BLOB,
                nat               BLOB,
                init_lattice      BLOB,
                init_frac_coords  BLOB,
                init_spg_num      INTEGER,
                init_spg_sym      TEXT,
                opt_lattice       BLOB,
                opt_frac_coords   BLOB,
                opt_spg_num       INTEGER,
                opt_spg_sym       TEXT,
                energy            REAL
            )
            """
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS ea_origin (
                id                    INTEGER PRIMARY KEY,
                generation            INTEGER NOT NULL,
                operation             TEXT NOT NULL,
                parent1_id            INTEGER,
                parent2_id            INTEGER,
                parent_attempt_count  INTEGER NOT NULL DEFAULT 0,
                FOREIGN KEY (id) REFERENCES records(id)
            )
            """
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS ea_info (
                generation      INTEGER PRIMARY KEY,
                population      INTEGER NOT NULL,
                n_crossover     INTEGER NOT NULL,
                n_permutation   INTEGER NOT NULL,
                n_strain        INTEGER NOT NULL,
                n_random        INTEGER NOT NULL,
                n_elite         INTEGER NOT NULL,
                crs_lat         TEXT NOT NULL,
                slct_func       TEXT NOT NULL,
                n_addition      INTEGER,
                n_elimination   INTEGER,
                n_substitution  INTEGER,
                min_comp        BLOB,
                max_comp        BLOB
            )
            """
        )
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS ea_elite (
                generation  INTEGER NOT NULL,
                id          INTEGER NOT NULL,
                PRIMARY KEY (generation, id),
                FOREIGN KEY (id) REFERENCES records(id)
            )
            """
        )
        conn.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_records_generation_id
            ON records(generation, id)
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