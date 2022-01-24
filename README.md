![cryspy_logo](./cryspy_fix-03.png)

# CrySPY
CrySPY is a crystal structure prediction tool written in Python.  
Document site is moved to https://tomoki-yamashita.github.io/CrySPY_doc/

## Latest version
version 0.10.2 (2022 January 24)

## Important changes
* version 0.10.2 (2022 January 24)
    - Added nrot: maximum number of times to rotate molecules in mol_bs mode
* version 0.10.1 (2021 September 30)
    - Bug fixed for numpy.random.seed in multiprocessing
* version 0.10.0 (2021 July 25)
    - Support PyXtal 0.2.9 or later
    - Step data and LAQA for QE
    - Upper and lower limits of energy for EA and BO
* version 0.9.2 (2021 Mar 18)
    - Support pymatgen major change (v2022)
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
- Python 3.8.x, 3.9.x (3.7.x may work)
- [COMBO](https://github.com/Tomoki-YAMASHITA/combo3 "COMBO")
- [PyXtal >= 0.2.2](https://pyxtal.readthedocs.io/en/latest "PyXtal")
- (PyXtal requires pymatgen) [pymatgen >= 2022.0.4](http://pymatgen.org "pymatgen")

See [CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc/installation/requirements/ "CrySPY document") in detail.

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
### CrySPY (software)
* T. Yamashita, S. Kanehira, N. Sato, H. Kino, H. Sawahata, T. Sato, F. Utsuno, K. Tsuda, T. Miyake, and T. Oguchi, Sci. Technol. Adv. Mater.:Methods **1**, 87 (2021).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2021.1943171


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
