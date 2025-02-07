from logging import getLogger
import warnings
from typing import Callable, Optional, Tuple

from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS, LBFGS, FIRE
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.structure import Structure


logger = getLogger('cryspy')


def opt_ase(
        work_path: str,
        struc: Structure,
        calculator: Callable,
        optimizer: str,
        symmetry: bool = True,
        fmax: float = 0.01,
        steps: int = 2000,
    ) -> Tuple[Optional[Structure], float, bool]:
    """
    Optimize the structure using ASE.

    Args:
        work_path (str): Path to the working directory.
        struc (Structure): Initial structure to optimize.
        calculator (Callable): Calculator to use for optimization.
        optimizer (str): Optimizer to use ('BFGS', 'LBFGS', 'FIRE').
        symmetry (bool, optional): Whether to use symmetry constraints. Default is True.
        fmax (float, optional): Maximum force. Default is 0.01.
        steps (int, optional): Number of steps. Default is 2000.

    Returns:
        Tuple[Optional[Structure], float, bool]:
            - Optimized structure
            - Total energy (eV/cell)
            - Whether the optimization converged
    """
    # ---------- atoms
    atoms = AseAtomsAdaptor.get_atoms(struc)

    # ---------- set calculator
    atoms.calc = calculator

    # ---------- set symmetry constraint
    if symmetry:
        atoms.set_constraint([FixSymmetry(atoms)])

    # ---------- set Filter
    cell_filter = FrechetCellFilter(atoms)

    # ---------- optimize structure
    oc_dict = {"BFGS": BFGS, "LBFGS": LBFGS, "FIRE": FIRE}
    oc = oc_dict[optimizer]
    opt = oc(cell_filter, trajectory=work_path+'opt.traj')
    with warnings.catch_warnings():
        # to avoid UserWarning: Only FixAtoms ~~ from Pymatgen
        warnings.simplefilter('ignore')
        converged = opt.run(fmax, steps)

    # ---------- result
    try:
        lattice = cell_filter.atoms.cell[:]
        species = cell_filter.atoms.get_chemical_symbols()
        coords = cell_filter.atoms.get_scaled_positions()
        opt_struc = Structure(lattice=lattice, species=species, coords=coords)
        energy = cell_filter.atoms.get_total_energy()    # eV/cell
    except Exception as e:
        opt_struc = None
        energy = np.nan
        logger.error(f'Optimization failed in {work_path}: {e}')

    # ---------- return
    return opt_struc, energy, converged