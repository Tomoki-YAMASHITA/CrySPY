import warnings

from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.optimize import BFGS, LBFGS, FIRE
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor


def opt_ase(
        work_path,
        struc,
        calculator,
        optimizer,
        symmetry=True,
        mask=None,
        fmax=0.01,
        steps=2000,
    ):
    # ---------- atoms
    atoms = AseAtomsAdaptor.get_atoms(struc)

    # ---------- set calculator
    atoms.calc = calculator

    # ---------- set symmetry constraint
    if symmetry:
        atoms.set_constraint([FixSymmetry(atoms)])

    # ---------- set Filter
    cell_filter = FrechetCellFilter(atoms, hydrostatic_strain=False, mask=mask)

    # ---------- optimize structure
    oc_dict = {"BFGS": BFGS, "LBFGS": LBFGS, "FIRE": FIRE}
    oc = oc_dict[optimizer]
    opt = oc(cell_filter, trajectory=work_path+'opt.traj')
    with warnings.catch_warnings():
        # to avoid UserWarning: Only FixAtoms ~~ from Pymatgen
        warnings.simplefilter('ignore')
        opt.run(fmax, steps)

    # ---------- result
    try:
        with warnings.catch_warnings():
            # to avoid UserWarning: Only FixAtoms ~~ from Pymatgen
            warnings.simplefilter('ignore')
            opt_struc = AseAtomsAdaptor.get_structure(cell_filter.atoms.copy())
            # ".copy()" <-- to prevent "TypeError: cannot pickle 'kimpy.model.PyModel' object" when using KIM
        energy = cell_filter.atoms.get_total_energy()    # eV/cell
    except Exception:
        opt_struc = None
        energy = np.nan

    # ---------- return
    return opt_struc, energy