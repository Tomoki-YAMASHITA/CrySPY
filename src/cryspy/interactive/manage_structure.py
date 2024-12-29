from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor

from ..IO import pkl_data


def get_ase_atoms(status: str, cid: int) -> Atoms:
    """
    Get ASE Atoms object based on the given status and structure ID.

    Args:
        status (str): Specify the status, either 'init' or 'opt'.
        cid (int): Specify the structure ID.

    Returns:
        Atoms: ASE Atoms object.

    Raises:
        ValueError: If the status is not 'init' or 'opt'.
        KeyError: If the specified structure ID is not found in the data.
    """
    # ---------- load structure data
    if status == 'init':
        struc_data = pkl_data.load_init_struc()
    elif status == 'opt':
        struc_data = pkl_data.load_opt_struc()
    else:
        raise ValueError('status should be init or opt')

    # ---------- get atoms
    struc = struc_data[cid]
    atoms = AseAtomsAdaptor.get_atoms(struc)

    # ---------- return
    return atoms
