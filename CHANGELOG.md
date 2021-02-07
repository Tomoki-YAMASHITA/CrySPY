# Changelog

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
