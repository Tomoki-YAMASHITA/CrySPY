![cryspy_logo](https://github.com/Tomoki-YAMASHITA/CrySPY/blob/master/logo/cryspy_fix-03.png)

[![PyPI version](https://badge.fury.io/py/csp-cryspy.svg)](https://badge.fury.io/py/csp-cryspy)
[![Downloads](https://static.pepy.tech/badge/csp-cryspy)](https://pepy.tech/project/csp-cryspy)

# CrySPY
CrySPY (pronounced as crispy) is a crystal structure prediction tool written in Python.  
Document: https://tomoki-yamashita.github.io/CrySPY_doc  
Questions and comments: https://github.com/Tomoki-YAMASHITA/CrySPY/discussions

## Latest version
version 1.3.0 (2024 May 31)

## News
- [2024 May 31] CrySPY 1.3.0 released.
    + There are important changes. See [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info)
- [2024 May 10] CrySPY 1.2.5 released.
    + bug fix for order_ef in out_results.py
- [2024 May 7] CrySPY 1.2.4 released.
    + bug fix
- [2023 October 21] CrySPY 1.2.3 released.
    + bug fix for MPI
- [2023 October 18] CrySPY 1.2.2 released.
    + [Enthalpy](https://tomoki-yamashita.github.io/CrySPY_doc/features/enthalpy/index.html)
- [2023 September 27] CrySPY 1.2.1 released.
    + bug fix for ASE interface
- [2023 July 10] CrySPY 1.2.0 released. Version information/version 1.2.0
    + Interface for ASE
    + Adoption of logging
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)
- [2023 June 14] CrySPY 1.1.1 released
    + bug fix
- [2023 May 16] CrySPY 1.1.0 released
    + MPI parallelization (optional)
    + New score of LAQA
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)
- [2023 March 16] CrySPY 1.0.0 released
    + CrySPY is available in PyPI, so you can install by pip (project name is csp-cryspy).
    + See also [Version information](https://tomoki-yamashita.github.io/CrySPY_doc/version_info) and [CHANGELOG](./CHANGELOG.md)


## System requirements
### Python
- Python >= 3.8
- [PyXtal >= 0.5.3](https://pyxtal.readthedocs.io/en/latest "PyXtal")

(optional)
- [PHYSBO](https://www.pasums.issp.u-tokyo.ac.jp/physbo/en/about "PHYSBO") (required if algo is BO)
- [DScribe](https://singroup.github.io/dscribe/latest/ "DScribe") (required if algo is BO)
- [mpi4py](https://mpi4py.readthedocs.io/en/stable "mpi4py")


See [CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc/installation/requirements/ "CrySPY document") in detail.

### Structure optimizer
At least one optimizer is required.

- [VASP](https://www.vasp.at "VASP") (tested with version 5.4.4)
- [QUANTUM ESPRESSO](http://www.quantum-espresso.org "Quantum ESPRESSO") (tested with version 6.x, version 5.x does not work)
- [OpenMX](http://www.openmx-square.org "OpenMX")
- [soiap](https://github.com/nbsato/soiap "soiap") (tested with version 0.2.2)
- [LAMMPS](http://lammps.sandia.gov "LAMMPS")
- [ASE](https://wiki.fysik.dtu.dk/ase "ASE")


## Document (English/Japanese)
[CrySPY document](https://tomoki-yamashita.github.io/CrySPY_doc "CrySPY documment")

## CrySPY Utility
[CrySPY Utility](https://github.com/Tomoki-YAMASHITA/CrySPY_utility "CrySPY Utility")

## Reference
### CrySPY (software)
* T. Yamashita, S. Kanehira, N. Sato, H. Kino, H. Sawahata, T. Sato, F. Utsuno, K. Tsuda, T. Miyake, and T. Oguchi, Sci. Technol. Adv. Mater. Meth. **1**, 87 (2021).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2021.1943171


### Bayesian optimization
* T. Yamashita, N. Sato, H. Kino, T. Miyake, K. Tsuda, and T. Oguchi, Phys. Rev. Mater. **2**, 013803 (2018).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013803

* N. Sato, T. Yamashita, T. Oguchi, K. Hukushima, and T. Miyake, Phys. Rev. Mater. **4**, 033801 (2020).
    - https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.4.033801

### Baysian optimization and evolutionary algorithm
* T. Yamashita, H. Kino, K. Tsuda, T. Miyake, and T. Oguchi, Sci. Technol. Adv. Mater. Meth. **2**, 67 (2022).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2022.2055987

### LAQA
* K.Terayama, T. Yamashita, T. Oguchi, and K. Tsuda, npj Comput. Mater. **4**, 32 (2018).
    - https://www.nature.com/articles/s41524-018-0090-y

* T. Yamashita and H. Sekine, Sci. Technol. Adv. Mater. Meth. **2**, 84 (2022).
    - https://www.tandfonline.com/doi/full/10.1080/27660400.2022.2059335


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team
