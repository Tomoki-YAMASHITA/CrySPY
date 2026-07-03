from enum import IntEnum
import sqlite3

import numpy as np
from ase import Atoms
from pymatgen.core import Structure

from ...util.struc_util import get_nat
from .convert import (
    array_to_blob,
    atoms_to_raw,
    blob_to_array,
    raw_to_struc,
    struc_to_raw,
)


class Status(IntEnum):
    WAITING = 0
    RUNNING = 1
    DONE = 2
    ERROR = 3
    NOT_CONVERGED = 4
    MINDIST = 5


def insert_init_struc(
    conn: sqlite3.Connection,
    record_id: int,
    struc: Structure,
    atype: tuple[str, ...],
    symprec: float,
    status: Status = Status.WAITING,
) -> int:
    """Insert an initial structure into the records table."""

    # ---------- space group
    try:
        init_spg_sym, init_spg_num = struc.get_space_group_info(
            symprec=symprec
        )
    except TypeError:
        init_spg_num = 0
        init_spg_sym = None

    # ---------- structure data
    raw_struc_dict = struc_to_raw(struc)
    atomic_numbers = raw_struc_dict["atomic_numbers"]
    lattice = raw_struc_dict["lattice"]
    frac_coords = raw_struc_dict["frac_coords"]
    nat = np.asarray(
        get_nat(struc, atype),
        dtype=np.uint16,
    )

    # ---------- insert into records table
    cursor = conn.execute(
        """
        INSERT INTO records (
            id,
            status,
            atomic_numbers,
            nat,
            init_lattice,
            init_frac_coords,
            init_spg_num,
            init_spg_sym,
            opt_lattice,
            opt_frac_coords,
            opt_spg_num,
            opt_spg_sym,
            energy
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            record_id,
            status,
            array_to_blob(atomic_numbers),
            array_to_blob(nat),
            array_to_blob(lattice),
            array_to_blob(frac_coords),
            init_spg_num,
            init_spg_sym,
            None,
            None,
            0,
            None,
            None,
        ),
    )

    # ---------- return
    return cursor.lastrowid


def select_record_ids(
    conn: sqlite3.Connection,
) -> list[int]:
    """Select all record IDs."""

    # ---------- select IDs
    cursor = conn.execute(
        """
        SELECT id
        FROM records
        ORDER BY id ASC
        """
    )

    # ---------- return
    return [record_id for record_id, in cursor]


def select_ids_by_status(
    conn: sqlite3.Connection,
) -> dict[Status, list[int]]:
    """Select record IDs grouped by status."""

    # ---------- initialize IDs
    ids_by_status = {status: [] for status in Status}

    # ---------- select IDs
    cursor = conn.execute(
        """
        SELECT id, status
        FROM records
        ORDER BY status ASC, id ASC
        """
    )
    for record_id, status in cursor:
        ids_by_status[Status(status)].append(record_id)

    # ---------- return
    return ids_by_status


def get_min_energy(
    conn: sqlite3.Connection,
) -> float | None:
    """Get the minimum energy."""

    # ---------- minimum energy
    cursor = conn.execute(
        """
        SELECT MIN(energy)
        FROM records
        WHERE status NOT IN (?, ?)
          AND energy IS NOT NULL
        """,
        (
            Status.WAITING,
            Status.RUNNING,
        ),
    )

    # ---------- return
    return cursor.fetchone()[0]


def count_results(
    conn: sqlite3.Connection,
    ids: list[int] = None,
    emin: float = None,
    emax: float = None,
) -> int:
    """Count calculated results."""

    # ---------- query
    query = """
        SELECT COUNT(*)
        FROM records
        WHERE status NOT IN (?, ?)
    """
    params = [
        Status.WAITING,
        Status.RUNNING,
    ]

    # ---------- IDs
    if ids is not None:
        placeholders = ', '.join('?' for _ in ids)
        query += f' AND id IN ({placeholders})'
        params.extend(ids)

    # ---------- energy range
    if emin is not None:
        query += ' AND energy >= ?'
        params.append(emin)
    if emax is not None:
        query += ' AND energy <= ?'
        params.append(emax)

    # ---------- return
    cursor = conn.execute(query, params)
    return cursor.fetchone()[0]


def select_results(
    conn: sqlite3.Connection,
    ids: list[int] = None,
    emin: float = None,
    emax: float = None,
    sort_by: str = 'energy',
    limit: int = None,
    tail: bool = False,
) -> list[tuple]:
    """Select calculated results."""

    # ---------- query
    query = """
        SELECT
            id,
            status,
            energy,
            init_spg_num,
            init_spg_sym,
            opt_spg_num,
            opt_spg_sym,
            nat
        FROM records
        WHERE status NOT IN (?, ?)
    """
    params = [
        Status.WAITING,
        Status.RUNNING,
    ]

    # ---------- IDs
    if ids is not None:
        placeholders = ', '.join('?' for _ in ids)
        query += f' AND id IN ({placeholders})'
        params.extend(ids)

    # ---------- energy range
    if emin is not None:
        query += ' AND energy >= ?'
        params.append(emin)
    if emax is not None:
        query += ' AND energy <= ?'
        params.append(emax)

    # ---------- order
    if sort_by == 'energy':
        has_energy_range = emin is not None or emax is not None
        if tail:
            if has_energy_range:
                query += ' ORDER BY energy DESC, id DESC'
            else:
                query += (
                    ' ORDER BY energy IS NULL DESC, '
                    'energy DESC, id DESC'
                )
        else:
            if has_energy_range:
                query += ' ORDER BY energy ASC, id ASC'
            else:
                query += (
                    ' ORDER BY energy IS NULL, '
                    'energy ASC, id ASC'
                )
    elif sort_by == 'id':
        if tail:
            query += ' ORDER BY id DESC'
        else:
            query += ' ORDER BY id ASC'
    else:
        raise ValueError(f'Unknown sort key: {sort_by}')

    # ---------- limit
    if limit is not None:
        query += ' LIMIT ?'
        params.append(limit)

    # ---------- select
    cursor = conn.execute(query, params)
    rows = cursor.fetchall()

    # ---------- restore ascending order for tail
    if tail:
        rows.reverse()

    # ---------- return
    return [
        (
            record_id,
            Status(status),
            energy,
            init_spg_num,
            init_spg_sym,
            opt_spg_num,
            opt_spg_sym,
            tuple(blob_to_array(nat).tolist()),
        )
        for (
            record_id,
            status,
            energy,
            init_spg_num,
            init_spg_sym,
            opt_spg_num,
            opt_spg_sym,
            nat,
        ) in rows
    ]


def select_struc(
    conn: sqlite3.Connection,
    record_id: int,
    initial: bool = False,
) -> dict | None:
    """Select structure data by ID."""

    # ---------- select structure
    cursor = conn.execute(
        """
        SELECT
            atomic_numbers,
            init_lattice,
            init_frac_coords,
            opt_lattice,
            opt_frac_coords
        FROM records
        WHERE id = ?
        """,
        (record_id,),
    )
    row = cursor.fetchone()

    # ---------- ID not found
    if row is None:
        return None

    # ---------- select initial or optimized structure
    if initial:
        lattice_blob = row[1]
        frac_coords_blob = row[2]
    else:
        lattice_blob = row[3]
        frac_coords_blob = row[4]

    # ---------- return
    return {
        'atomic_numbers': blob_to_array(row[0]),
        'lattice': (
            None
            if lattice_blob is None
            else blob_to_array(lattice_blob)
        ),
        'frac_coords': (
            None
            if frac_coords_blob is None
            else blob_to_array(frac_coords_blob)
        ),
    }


def update_opt_struc(
    conn: sqlite3.Connection,
    record_id: int,
    opt_atoms: Atoms,
    energy: float,
    symprec: float,
) -> None:
    """Update an optimized structure in the records table."""

    # ---------- structure data
    raw_struc_dict = atoms_to_raw(opt_atoms)
    lattice = raw_struc_dict["lattice"]
    frac_coords = raw_struc_dict["frac_coords"]

    # ---------- pymatgen Structure for space group
    opt_struc = raw_to_struc(raw_struc_dict)

    # ---------- space group
    try:
        opt_spg_sym, opt_spg_num = opt_struc.get_space_group_info(
            symprec=symprec
        )
    except TypeError:
        opt_spg_num = 0
        opt_spg_sym = None

    # ---------- update records table
    conn.execute(
        """
        UPDATE records
        SET
            opt_lattice = ?,
            opt_frac_coords = ?,
            opt_spg_num = ?,
            opt_spg_sym = ?,
            energy = ?
        WHERE id = ?
        """,
        (
            array_to_blob(lattice),
            array_to_blob(frac_coords),
            opt_spg_num,
            opt_spg_sym,
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