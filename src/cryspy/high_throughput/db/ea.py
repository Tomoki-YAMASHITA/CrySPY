import sqlite3


def insert_ea_origin(
    conn: sqlite3.Connection,
    record_id: int,
    generation: int,
    operation: str,
    parent1_id: int | None = None,
    parent2_id: int | None = None,
    parent_attempt_count: int = 0,
) -> None:
    """Insert EA origin information."""

    # ---------- insert ea_origin table
    conn.execute(
        """
        INSERT INTO ea_origin (
            id,
            generation,
            operation,
            parent1_id,
            parent2_id,
            parent_attempt_count
        )
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (
            record_id,
            generation,
            operation,
            parent1_id,
            parent2_id,
            parent_attempt_count,
        ),
    )


def insert_ea_info(
    conn: sqlite3.Connection,
    generation: int,
    population: int,
    n_crossover: int,
    n_permutation: int,
    n_strain: int,
    n_random: int,
    n_elite: int,
    crs_lat: str,
    slct_func: str,
    n_addition: int | None = None,
    n_elimination: int | None = None,
    n_substitution: int | None = None,
    min_comp: bytes | None = None,
    max_comp: bytes | None = None,
) -> None:
    """Insert EA generation information."""

    # ---------- insert ea_info table
    conn.execute(
        """
        INSERT INTO ea_info (
            generation,
            population,
            n_crossover,
            n_permutation,
            n_strain,
            n_random,
            n_elite,
            crs_lat,
            slct_func,
            n_addition,
            n_elimination,
            n_substitution,
            min_comp,
            max_comp
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            generation,
            population,
            n_crossover,
            n_permutation,
            n_strain,
            n_random,
            n_elite,
            crs_lat,
            slct_func,
            n_addition,
            n_elimination,
            n_substitution,
            min_comp,
            max_comp,
        ),
    )


def insert_ea_elites(
    conn: sqlite3.Connection,
    generation: int,
    elite_ids: list[int],
) -> None:
    """Insert EA elite information."""

    # ---------- insert ea_elite table
    conn.executemany(
        """
        INSERT INTO ea_elite (
            generation,
            id
        )
        VALUES (?, ?)
        """,
        [
            (
                generation,
                elite_id,
            )
            for elite_id in elite_ids
        ],
    )