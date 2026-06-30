from enum import IntEnum
import sqlite3

import numpy as np
from ase import Atoms
from pymatgen.core import Structure

from .convert import array_to_blob, atoms_to_raw, blob_to_array, struc_to_raw


class Status(IntEnum):
    WAITING = 0
    RUNNING = 1
    DONE = 2
    ERROR = 3
    NOT_CONVERGED = 4
    MINDIST = 5


def insert_init_struc(
    conn: sqlite3.Connection,
    struc: Structure,
    nat: np.ndarray,
    status: Status = Status.WAITING,
) -> int:
    """Insert an initial structure into the records table."""

    # ---------- structure data
    raw_struc_dict = struc_to_raw(struc)
    atomic_numbers = raw_struc_dict["atomic_numbers"]
    lattice = raw_struc_dict["lattice"]
    frac_coords = raw_struc_dict["frac_coords"]
    nat = np.asarray(nat, dtype=np.uint16)

    # ---------- insert into records table
    cursor = conn.execute(
        """
        INSERT INTO records (
            status,
            atomic_numbers,
            nat,
            init_lattice,
            init_frac_coords,
            opt_lattice,
            opt_frac_coords,
            energy
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            status,
            array_to_blob(atomic_numbers),
            array_to_blob(nat),
            array_to_blob(lattice),
            array_to_blob(frac_coords),
            None,
            None,
            None,
        ),
    )

    # ---------- return
    return cursor.lastrowid


# select_init_strucはもう不要
# select_next_strucがある。

# def select_init_struc(
#     conn: sqlite3.Connection,
#     record_id: int,
# ) -> dict | None:
#     """Select initial structure data by id."""

#     # ---------- select from records table
#     cursor = conn.execute(
#         """
#         SELECT
#             atomic_numbers,
#             nat,
#             init_lattice,
#             init_frac_coords
#         FROM records
#         WHERE id = ?
#         """,
#         (record_id,),
#     )
#     row = cursor.fetchone()

#     # ---------- return
#     if row is None:
#         return None
#     return {
#         "atomic_numbers": blob_to_array(row[0]),
#         "nat": blob_to_array(row[1]),
#         "lattice": blob_to_array(row[2]),
#         "frac_coords": blob_to_array(row[3]),
#     }


def update_opt_struc(
    conn: sqlite3.Connection,
    record_id: int,
    opt_atoms: Atoms,
    energy: float,
) -> None:
    """Update an optimized structure in the records table."""

    # ---------- structure data
    raw_struc_dict = atoms_to_raw(opt_atoms)
    lattice = raw_struc_dict["lattice"]
    frac_coords = raw_struc_dict["frac_coords"]

    # ---------- update records table
    conn.execute(
        """
        UPDATE records
        SET
            opt_lattice = ?,
            opt_frac_coords = ?,
            energy = ?
        WHERE id = ?
        """,
        (
            array_to_blob(lattice),
            array_to_blob(frac_coords),
            energy,
            record_id,
        ),
    )


def claim_next_struc(
    conn: sqlite3.Connection,
) -> tuple[int, dict] | None:
    """Atomically claim the next structure to calculate."""

    # ---------- begin write transaction
    conn.execute("BEGIN IMMEDIATE")

    try:
        # ---------- select from records table
        cursor = conn.execute(
            """
            SELECT
                id,
                atomic_numbers,
                nat,
                init_lattice,
                init_frac_coords
            FROM records
            WHERE status = ?
            ORDER BY id ASC
            LIMIT 1
            """,
            (Status.WAITING,),
        )
        row = cursor.fetchone()

        # ---------- if no record found
        if row is None:
            conn.commit()
            return None

        # ---------- update status
        record_id = row[0]
        conn.execute(
            """
            UPDATE records
            SET status = ?
            WHERE id = ?
            """,
            (
                Status.RUNNING,
                record_id,
            ),
        )

        # ---------- structure data
        raw_struc_dict = {
            "atomic_numbers": blob_to_array(row[1]),
            "nat": blob_to_array(row[2]),
            "lattice": blob_to_array(row[3]),
            "frac_coords": blob_to_array(row[4]),
        }

        # ---------- commit
        conn.commit()

    except Exception:
        # ---------- rollback
        conn.rollback()
        raise

    # ---------- return
    return record_id, raw_struc_dict


def reset_running_struc(
    conn: sqlite3.Connection,
) -> int:
    """Reset running structures to waiting."""

    # ---------- update status
    cursor = conn.execute(
        """
        UPDATE records
        SET status = ?
        WHERE status = ?
        """,
        (
            Status.WAITING,
            Status.RUNNING,
        ),
    )

    # ---------- return
    return cursor.rowcount


def update_status(
    conn: sqlite3.Connection,
    record_id: int,
    status: Status,
) -> None:
    """Update status in the records table."""

    # ---------- update records table
    conn.execute(
        """
        UPDATE records
        SET status = ?
        WHERE id = ?
        """,
        (
            status,
            record_id,
        ),
    )