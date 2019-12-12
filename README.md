# CrySPY
CrySPY is a crystal structure prediction tool written in Python.

You can automatically do the following things with CrySPY.
* Structure generation
* Submitting jobs for structure optimization
* Collecting data for structure optimization
* Selecting candidates using machine learning

You can also use CrySPY partly for the purpose of only struture generation and descriptor calculation for crystal structures.

## Latest version
version 0.7.0 (2018 December 5)

## Important changes
* version 0.7.0
    - Evolutionary algorithm is now available
* version 0.6.2
    - LAMMPS can be used in CrySPY
* version 0.6.0
    - LAQA is now available
    - Changed the data format of init_struc_data and opt_struc_data from list to dict

## System requirements
### Python
- Python 2.7.x
- [COMBO](https://github.com/tsudalab/combo "COMBO")
- numpy
- pandas
- [pymatgen](http://pymatgen.org "pymatgen")

### Structure optimizer
At least either one of

- [VASP](https://www.vasp.at "VASP") (tested with version 5.4.1)
- [Quantum ESPRESSO](http://www.quantum-espresso.org "Quantum ESPRESSO") (tested with version 6.1, version 5.x does not work)
- [soiap](https://github.com/nbsato/soiap "soiap") (tested with version 0.2.2)
- [LAMMPS](http://lammps.sandia.gov "LAMMPS")

### Others
- [find_wy](https://github.com/nim-hrkn/find_wy "find_wy")

## Document
[CrySPY document](https://tomoki-yamashita.github.io/CrySPY "CrySPY documment")

## Tutorial (written in Japanese)
[チュートリアルと解説](https://tomoki-yamashita.github.io/cryspy/tutorial/outline.html "tutorial")

## Reference
### Bayesian optimization
T. Yamashita, N. Sato, H. Kino, T. Miyake, K. Tsuda, and T. Oguchi, Phys. Rev. Mater. **2**, 013803 (2018).  
https://link.aps.org/doi/10.1103/PhysRevMaterials.2.013803

### LAQA
K.Terayama, T. Yamashita, T. Oguchi, and K. Tsuda, npj Comput. Mater. **4**, 32 (2018).  
https://www.nature.com/articles/s41524-018-0090-y


## License
CrySPY is distributed under the MIT License.  
Copyright (c) 2018 CrySPY Development Team
