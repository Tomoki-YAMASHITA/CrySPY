from io import BytesIO

from ase import Atoms
import numpy as np
from pymatgen.core import Lattice, Structure


def array_to_blob(arr: np.ndarray) -> bytes:
    """Convert a NumPy array to a SQLite BLOB (.npy binary)."""
    bio = BytesIO()
    np.save(bio, arr, allow_pickle=False)
    return bio.getvalue()


def blob_to_array(blob: bytes) -> np.ndarray:
    """Convert a SQLite BLOB (.npy binary) to a NumPy array."""
    bio = BytesIO(blob)
    return np.load(bio, allow_pickle=False)


def raw_to_struc(
    raw_struc_dict: dict,
) -> Structure:
    """Convert a raw structure dictionary to a pymatgen Structure."""

    return Structure(
        lattice=Lattice(raw_struc_dict["lattice"]),
        species=raw_struc_dict["atomic_numbers"],
        coords=raw_struc_dict["frac_coords"],
        coords_are_cartesian=False,
    )


def struc_to_raw(
    struc: Structure,
) -> dict:
    """Convert a pymatgen Structure to a raw structure dictionary."""

    return {
        "atomic_numbers": np.array(
            [site.specie.Z for site in struc],
            dtype=np.uint8,
        ),
        "lattice": np.array(
            struc.lattice.matrix,
            dtype=np.float32,
        ),
        "frac_coords": np.array(
            struc.frac_coords,
            dtype=np.float32,
        ),
    }


def raw_to_atoms(
    raw_struc_dict: dict,
) -> Atoms:
    """Convert a raw structure dictionary to ASE Atoms."""

    return Atoms(
        numbers=raw_struc_dict["atomic_numbers"],
        cell=raw_struc_dict["lattice"],
        scaled_positions=raw_struc_dict["frac_coords"],
        pbc=True,
    )


def atoms_to_raw(
    atoms: Atoms,
) -> dict:
    """Convert ASE Atoms to a raw structure dictionary."""

    return {
        "atomic_numbers": np.array(
            atoms.get_atomic_numbers(),
            dtype=np.uint8,
        ),
        "lattice": np.array(
            atoms.cell.array,
            dtype=np.float32,
        ),
        "frac_coords": np.array(
            atoms.get_scaled_positions(),
            dtype=np.float32,
        ),
    }