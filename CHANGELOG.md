# CHANGELOG
## [1.3.0] - 2024-5-31
### Important change
#### Common
- working directory name  
  work000000 --> work0
- We used to pickle data by grouping several data into tuples, but we changed it to pickle each item individually.  
  For example, rs_id_data.pkl --> id_queueing.pkl and id_running.pkl
#### BO
- Using cal_fingerprint program is obsolete. [dscribe](https://singroup.github.io/dscribe/latest) is required instead.
  + `fppath` and `fp_rmin` in cryspy.in are obsolete.
- Changed the Bayesian optimization library from COMBO to [PHYSBO](https://www.pasums.issp.u-tokyo.ac.jp/physbo/en/about)
### Fixed
#### soiap
- support for recent pymatgen
### Added
- Random structure generation and structure generation by EA are now available as libraries.
### for developer
- We stopped using global variables (rin), now uses dataclass for input data.
- Many of the input variables were lists, but we changed them to tuples.

## [1.2.5] - 2024-5-10
### Fixed
- bug fix for order_ef in out_results.py

## [1.2.4] - 2024-5-7
### Fixed
- bug fix
### Changed
- default value of cls_lat: equal –> random
### Added
- test version of variable composition EA (EA-vc). only binary system for now.

## [1.2.3] - 2023-10-21
### Fixed
- bug fix for MPI

## [1.2.2] - 2023-10-18
### Added
- enthalpy

## [1.2.1] - 2023-09-27
### Fixed
- bug fix for ASE interface

## [1.2.0] - 2023-07-10
### Changed
- adopts logging for output
### Added
- ASE interface

## [1.1.1] - 2023-06-14
### Fixed
- bug fix for gen_struc: delete spg_error

## [1.1.0] - 2023-05-16
### Changed
- backup (keep history)
### Added
- MPI parallelization
- New score of LAQA

## [1.0.0] - 2023-03-16
### Changed
- standard output --> log_cryspy
- standard error --> err_cryspy
- cryspy.out is obsoleted
- COMBO --> optional
- mindist
  + mindist can be omitted in cryspy.in
  + mindist_ea is obsoleted
- moved examples to CrySPY utility
- moved cal_fingerprint program to CrySPY utility
- directory tree (gen_struc, utility, RS, EA)
- CrySPY stops once before going to next selection or next generation

### Added
- excecutable script, cryspy
- auto and manual backup
- mindist_mol_bs and mindist_mol_bs_factor in cryspy.in
- fwpath and fppath in cryspy.in
- new calc_mode: ext


## [0.10.3] - 2022-05-17
### Changed
- update LAMMPS example
- default vol_factor value: 1.1

### Fixed
- bug fix for LAMMPS IO
- parse force data in QE

## [0.10.2] - 2022-01-24
### Added
- nrot: maximum number of times to rotate molecules in mol_bs

## [0.10.1] - 2021-09-30
### Fixed
- bug fix for np.random.seed in multiprocessing

## [0.10.0] - 2021-07-25
### Added
- step_data for QE
- LAQA for QE
- Upper and lower limits of energy in EA and BO

### Changed
- support PyXtal 0.2.9 or later

### Fixed
- bug fix for recalc

## [0.9.2] - 2021-03-18
### Changed
- support pymatgen major change (v2022)

## [0.9.1] - 2021-02-25
### Fixed
- minor bug fix for bo_status

## [0.9.0] - 2021-02-07
### Added
- Interface with OpenMX
- Employ PyXtal
  + Molecular crystal structure generation
- Scale volume
- LAQA with soiap

### Changed
- find_wy is optional
- In cryspy.in, [lattice] section --> [structure] section
  + natot: [basic] –> [structure]
  + atype: [basic] –> [structure]
  + nat: [basic] –> [structure]
  + maxcnt: [option] –> [structure]
  + symprec: [option] –> [structure]
  + spgnum: [option] –> [structure]

## [0.8.0] - 2020-02-16
### Added
- CrySPY logo
- recaluculation
- manual select in BO
### Changed
- Migrated to Python 3
- Several variable names
- Several data formats
- Unit of energy in output: eV/cell --> eV/atom
- No. of working directories corresponds to structure ID

## [0.7.0] - 2018-12-05
### Added
- Evolutionary algorithm is now available
- Added lock_cryspy system
### Fixed
- Fixed minor bugs

## [0.6.4] - 2018-08-20
### Fixed
- Fixed a bug for spgnum in cryspy.in/read_input.py

## [0.6.3] - 2018-07-02
### Added
- Added LAMMPS example
### Changed
- In generating structures (including appending structures), CrySPY always stops before submitting jobs.

## [0.6.2] - 2018-06-26
### Added
- LAMMPS can be used in CrySPY

## [0.6.1] - 2018-03-01
### Fixed
- Fixed a bug for job control in RS

## [0.6.0] - 2018-02-20
### Added
- LAQA is now available

### Changed
- Changed the data format of init_struc_data and opt_struc_data from list to dict

### Fixed
- Fixed a minor bug for BO
