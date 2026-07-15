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


def select_ea_elite_ids(
    conn: sqlite3.Connection,
    generation: int,
) -> list[int]:
    """Select EA elite IDs."""

    # ---------- select elite IDs
    cursor = conn.execute(
        """
        SELECT id
        FROM ea_elite
        WHERE generation = ?
        ORDER BY id ASC
        """,
        (generation,),
    )

    # ---------- return
    return [elite_id for (elite_id,) in cursor]


def select_latest_ea_generation(
    conn: sqlite3.Connection,
) -> int | None:
    """Select latest EA generation."""

    # ---------- select latest generation
    cursor = conn.execute(
        """
        SELECT MAX(generation)
        FROM ea_info
        """
    )

    # ---------- return
    return cursor.fetchone()[0]


def select_ea_info(
    conn: sqlite3.Connection,
    generation: int,
) -> tuple | None:
    """Select EA generation information."""

    # ---------- select ea_info
    cursor = conn.execute(
        """
        SELECT
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
        FROM ea_info
        WHERE generation = ?
        """,
        (generation,),
    )

    # ---------- return
    return cursor.fetchone()


def select_ea_origins(
    conn: sqlite3.Connection,
    generation: int,
) -> list[tuple]:
    """Select EA origin information."""

    # ---------- select ea_origin
    cursor = conn.execute(
        """
        SELECT
            id,
            generation,
            operation,
            parent1_id,
            parent2_id,
            parent_attempt_count
        FROM ea_origin
        WHERE generation = ?
        ORDER BY id ASC
        """,
        (generation,),
    )

    # ---------- return
    return cursor.fetchall()


def select_ea_generations(
    conn: sqlite3.Connection,
    min_generation: int | None = None,
    max_generation: int | None = None,
) -> list[int]:
    """Select EA generations."""

    # ---------- query
    query = """
        SELECT generation
        FROM ea_info
        WHERE 1 = 1
    """
    params = []

    # ---------- generation range
    if min_generation is not None:
        query += ' AND generation >= ?'
        params.append(min_generation)
    if max_generation is not None:
        query += ' AND generation <= ?'
        params.append(max_generation)

    # ---------- order
    query += ' ORDER BY generation ASC'

    # ---------- select generations
    cursor = conn.execute(
        query,
        params,
    )

    # ---------- return
    return [
        generation
        for (generation,) in cursor
    ]


def delete_ea_origins_after_generation(
    conn: sqlite3.Connection,
    generation: int,
) -> None:
    """Delete EA origins after generation."""

    # ---------- delete ea_origin
    conn.execute(
        """
        DELETE FROM ea_origin
        WHERE generation > ?
        """,
        (generation,),
    )


def delete_ea_elites_after_generation(
    conn: sqlite3.Connection,
    generation: int,
) -> None:
    """Delete EA elites after generation."""

    # ---------- delete ea_elite
    conn.execute(
        """
        DELETE FROM ea_elite
        WHERE generation > ?
        """,
        (generation,),
    )