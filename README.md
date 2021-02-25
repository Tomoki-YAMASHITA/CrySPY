![cryspy_logo](./cryspy_fix-03.png)

# CrySPY
CrySPY is a crystal structure prediction tool written in Python.  
Document site is moved to https://tomoki-yamashita.github.io/CrySPY_doc/

## Latest version
version 0.9.1 (2021 Feb 25)

## Important changes
* version 0.9.0 (2021 Feb 7)
    - Interfaced with OpenMX
    - Employ PyXtal library to generate initial structures
    - If you use PyXtal (default), find_wy program is not required
    - LAQA can be used with soiap
    - Change the section name in the input file: [lattice] section –> [structure] section
    - Several input variables move to [structure] section
        + natot: [basic] –> [structure]
        + atype: [basic] –> [structure]
        + nat: [basic] –> [structure]
        + maxcnt: [option] –> [structure]
        + symprec: [option] –> [structure]
        + spgnum: [option] –> [structure]
    - New features
        + Molecular crystal structure generation
        + Scale volume
* version 0.8.0
    - Migrated to Python 3
    - Several variable names
    - Several data formats
    - Unit of energy in output: eV/cell --> eV/atom
    - No. of working directories corresponds to structure ID
* version 0.7.0
    - Evolutionary algorithm is now available
* version 0.6.2
    - LAMMPS can be used in CrySPY
* version 0.6.0
    - LAQA is now available
    - Changed the data format of init_struc_data and opt_struc_data from list to dict

## System requirements
### Python
- Python 3.8.x (3.7.x may work, PyXtal (numba) does not support Python 3.9)
- [COMBO](https://github.com/tsudalab/combo3 "COMBO")
- [pymatgen](http://pymatgen.org "pymatgen")
- [PyXtal >= 0.1.6](https://pyxtal.readthedocs.io/en/latest "PyXtal")

PyXtal requires SciPy, but the latest version of SciPy (v1.6.0) might include a bug for deepcopy.
You should use SciPy v1.5.4 for a while (2021 Feb 7).

### Structure optimizer
At least one optimizer is required.

- [VASP](https://www.vasp.at "VASP") (tested with version 5.4.4)
- [QUANTUM ESPRESSO](http://www.quantum-espresso.org "Quantum ESPRESSO") (tested with version 6.x, version 5.x does not work)
- [OpenMX](http://www.openmx-square.org "OpenMX")
- [soiap](https://github.com/nbsato/soiap "soiap") (tested with version 0.2.2)
- [LAMMPS](http://lammps.sandia.gov "LAMMPS")

### Others
- [find_wy](https://github.com/nim-hrkn/find_wy "find_wy"): find_wy can randomly select a combination of Wyckoff positions for a given chemical composition and space group. (optional)

## Document (English/Japanese)
[CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc "CrySPY documment")

## Google group
[Google gruop of CrySPY](https://groups.google.com/forum/#!forum/cryspy-user "Google group")



## Reference
### Bayesian optimization
* T. Yamashita, N. Sato, H. Kino, T. Miyake, K. Tsuda, and T. Oguchi, Phys. Rev. Mater. **2**, 013803 (2018).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013803

* N. Sato, T. Yamashita, T. Oguchi, K. Hukushima, and T. Miyake, Phys. Rev. Mater. **4**, 033801 (2020).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.033801


### LAQA
* K.Terayama, T. Yamashita, T. Oguchi, and K. Tsuda, npj Comput. Mater. **4**, 32 (2018).
    - https://www.nature.com/articles/s41524-018-0090-y


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team
