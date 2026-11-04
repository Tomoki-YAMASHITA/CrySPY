import sqlite3

import numpy as np

from .convert_db import (
    array_to_blob,
    blob_to_array,
    json_to_phase_diagram,
    phase_diagram_to_json,
)


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
    crs_lat: str | None,
    slct_func: str,
    n_addition: int | None = None,
    n_elimination: int | None = None,
    n_substitution: int | None = None,
    min_comp: tuple[float, ...] | None = None,
    max_comp: tuple[float, ...] | None = None,
    ll_nat: tuple[int, ...] | None = None,
    ul_nat: tuple[int, ...] | None = None,
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
            max_comp,
            ll_nat,
            ul_nat
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            array_to_blob(np.asarray(min_comp, dtype=float))
            if min_comp is not None else None,
            array_to_blob(np.asarray(max_comp, dtype=float))
            if max_comp is not None else None,
            array_to_blob(np.asarray(ll_nat, dtype=int))
            if ll_nat is not None else None,
            array_to_blob(np.asarray(ul_nat, dtype=int))
            if ul_nat is not None else None,
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


def delete_ea_elites(
    conn: sqlite3.Connection,
    generation: int,
) -> None:
    """Delete elite information for the given generation."""

    # ---------- delete elite IDs
    conn.execute(
        """
        DELETE FROM ea_elite
        WHERE generation = ?
        """,
        (generation,),
    )


def save_ea_hdist(
    conn: sqlite3.Connection,
    generation: int,
    hdist: dict[int, float],
) -> None:
    """Save hull distances for an EA generation."""

    # ---------- delete existing hull distances
    conn.execute(
        """
        DELETE FROM ea_hdist
        WHERE generation = ?
        """,
        (generation,),
    )

    # ---------- insert hull distances
    conn.executemany(
        """
        INSERT INTO ea_hdist (
            generation,
            id,
            hdist
        )
        VALUES (?, ?, ?)
        """,
        [
            (
                generation,
                cid,
                distance,
            )
            for cid, distance in hdist.items()
        ],
    )


def save_ea_phase_diagram(
    conn: sqlite3.Connection,
    generation: int,
    phase_diagram,
) -> None:
    """Save a phase diagram for an EA generation."""

    # ---------- save phase diagram
    conn.execute(
        """
        INSERT INTO ea_phase_diagram (
            generation,
            phase_diagram
        )
        VALUES (?, ?)
        ON CONFLICT (generation) DO UPDATE SET
            phase_diagram = excluded.phase_diagram
        """,
        (
            generation,
            phase_diagram_to_json(phase_diagram),
        ),
    )


def delete_ea_hull_data(
    conn: sqlite3.Connection,
    generation: int,
) -> None:
    """Delete convex hull information for the given generation."""

    # ---------- delete hull distances
    conn.execute(
        """
        DELETE FROM ea_hdist
        WHERE generation = ?
        """,
        (generation,),
    )

    # ---------- delete phase diagram
    conn.execute(
        """
        DELETE FROM ea_phase_diagram
        WHERE generation = ?
        """,
        (generation,),
    )

    # ---------- invalidate latest parent elites
    # Older generations retain their recorded elite history.
    if generation == select_latest_ea_generation(conn):
        delete_ea_elites(conn, generation)


def has_ea_phase_diagram(
    conn: sqlite3.Connection,
    generation: int,
) -> bool:
    """Return whether a phase diagram is saved for the generation."""

    # ---------- check existence
    row = conn.execute(
        """
        SELECT 1
        FROM ea_phase_diagram
        WHERE generation = ?
        """,
        (generation,),
    ).fetchone()
    return row is not None


def select_ea_phase_diagram(
    conn: sqlite3.Connection,
    generation: int,
):
    """Select a phase diagram for an EA generation."""

    # ---------- select phase diagram
    cursor = conn.execute(
        """
        SELECT phase_diagram
        FROM ea_phase_diagram
        WHERE generation = ?
        """,
        (generation,),
    )
    row = cursor.fetchone()

    # ---------- no phase diagram
    if row is None:
        return None

    # ---------- return
    return json_to_phase_diagram(row[0])


def select_ea_hdist(
    conn: sqlite3.Connection,
    generation: int,
) -> dict[int, float]:
    """Select hull distances for an EA generation."""

    # ---------- select hull distances
    cursor = conn.execute(
        """
        SELECT id, hdist
        FROM ea_hdist
        WHERE generation = ?
        ORDER BY id ASC
        """,
        (generation,),
    )

    # ---------- return
    return dict(cursor)


def select_ea_hdist_rows(
    conn: sqlite3.Connection,
    generation: int,
) -> list[tuple[int, int, int, float]]:
    """Select hull distances with evaluation and structure generations."""

    # ---------- select hull distances and structure generations
    cursor = conn.execute(
        """
        SELECT h.generation, r.generation, h.id, h.hdist
        FROM ea_hdist AS h
        JOIN records AS r ON r.id = h.id
        WHERE h.generation = ?
        ORDER BY h.id ASC
        """,
        (generation,),
    )

    # ---------- return
    return cursor.fetchall()


def select_ea_elite_rows(
    conn: sqlite3.Connection,
    generation: int,
) -> list[tuple[int, int, int]]:
    """Select elite IDs with elite and structure generations."""

    # ---------- select elite IDs and structure generations
    cursor = conn.execute(
        """
        SELECT e.generation, r.generation, e.id
        FROM ea_elite AS e
        JOIN records AS r ON r.id = e.id
        WHERE e.generation = ?
        ORDER BY e.id ASC
        """,
        (generation,),
    )

    # ---------- return
    return cursor.fetchall()


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


def select_ea_phase_diagram_status(
    conn: sqlite3.Connection,
    min_generation: int | None = None,
    max_generation: int | None = None,
    limit: int = 0,
    tail: bool = False,
) -> list[tuple[int, str, int]]:
    """Select PD storage status and hull distance counts without loading PD data."""

    # ---------- query
    query = """
        SELECT
            i.generation,
            CASE WHEN p.generation IS NULL THEN 'Missing' ELSE 'Saved' END,
            (
                SELECT COUNT(*)
                FROM ea_hdist AS h
                WHERE h.generation = i.generation
            )
        FROM ea_info AS i
        LEFT JOIN ea_phase_diagram AS p ON p.generation = i.generation
        WHERE 1 = 1
    """
    params = []

    # ---------- generation range
    if min_generation is not None:
        query += ' AND i.generation >= ?'
        params.append(min_generation)
    if max_generation is not None:
        query += ' AND i.generation <= ?'
        params.append(max_generation)

    # ---------- order
    if tail:
        query += ' ORDER BY i.generation DESC'
    else:
        query += ' ORDER BY i.generation ASC'

    # ---------- limit
    if limit > 0:
        query += ' LIMIT ?'
        params.append(limit)

    # ---------- select storage status
    rows = conn.execute(query, params).fetchall()

    # ---------- restore ascending order for tail
    if tail:
        rows.reverse()

    # ---------- return
    return rows


def select_latest_ea_phase_diagram_generation(
    conn: sqlite3.Connection,
) -> int | None:
    """Select latest generation with a saved phase diagram."""

    # ---------- select latest generation
    cursor = conn.execute(
        """
        SELECT MAX(generation)
        FROM ea_phase_diagram
        """
    )

    # ---------- return
    return cursor.fetchone()[0]


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
            max_comp,
            ll_nat,
            ul_nat
        FROM ea_info
        WHERE generation = ?
        """,
        (generation,),
    )

    # ---------- return
    return cursor.fetchone()


def select_ea_composition_bounds(
    conn: sqlite3.Connection,
    generation: int,
) -> tuple | None:
    """Select composition bounds for the given generation."""

    # ---------- select composition bounds
    row = conn.execute(
        """
        SELECT min_comp, max_comp
        FROM ea_info
        WHERE generation = ?
        """,
        (generation,),
    ).fetchone()
    if row is None:
        return None

    # ---------- return
    min_comp, max_comp = row
    return (
        tuple(blob_to_array(min_comp).tolist()) if min_comp is not None else None,
        tuple(blob_to_array(max_comp).tolist()) if max_comp is not None else None,
    )


def select_ea_parent_settings(
    conn: sqlite3.Connection,
    generation: int,
) -> tuple | None:
    """Select constraints and elite count for parent selection."""

    # ---------- select parent settings
    row = conn.execute(
        """
        SELECT ll_nat, ul_nat, n_elite, min_comp, max_comp
        FROM ea_info
        WHERE generation = ?
        """,
        (generation,),
    ).fetchone()
    if row is None:
        return None

    # ---------- return
    ll_nat, ul_nat, n_elite, min_comp, max_comp = row
    return (
        tuple(blob_to_array(ll_nat).tolist()) if ll_nat is not None else None,
        tuple(blob_to_array(ul_nat).tolist()) if ul_nat is not None else None,
        n_elite,
        tuple(blob_to_array(min_comp).tolist()) if min_comp is not None else None,
        tuple(blob_to_array(max_comp).tolist()) if max_comp is not None else None,
    )


def select_ea_algorithm(
    conn: sqlite3.Connection,
) -> str | None:
    """Select EA algorithm from initial generation information."""

    # ---------- select initial generation
    cursor = conn.execute(
        """
        SELECT n_addition, n_elimination, n_substitution
        FROM ea_info
        WHERE generation = 1
        """
    )
    row = cursor.fetchone()

    # ---------- no EA information
    if row is None:
        return None

    # ---------- algorithm
    return 'EA-vc' if any(value is not None for value in row) else 'EA'


def select_ea_origins(
    conn: sqlite3.Connection,
    generation: int,
) -> list[tuple]:
    """Select EA origin information."""

    # ---------- select ea_origin
    cursor = conn.execute(
        """
        SELECT
            generation,
            id,
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
