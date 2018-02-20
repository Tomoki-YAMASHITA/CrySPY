#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from pymatgen import Structure
from pymatgen.core.units import Length
from pymatgen.io.cif import CifWriter


def get_natot(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if 'nat' in line:
            nat = int(line.split()[-1])
            break

    return nat


def extract_cell_parameters(filename):
    #---------- last CELL_PARAMETERS
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_cell = None
    for i, line in enumerate(reversed(lines)):
        if 'CELL_PARAMETERS' in line:
            ibegin = len(lines) - 1 - i
            iend = ibegin + 4
            lines_cell = lines[ibegin:iend]
            break

    return lines_cell


def extract_atomic_positions(filename, natot):
    #---------- last ATOMIC_POSITIONS
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_atom = None
    for i, line in enumerate(reversed(lines)):
        if 'ATOMIC_POSITIONS' in line:
            ibegin = len(lines) - 1 - i
            iend = ibegin + natot + 1
            lines_atom = lines[ibegin:iend]
            break

    return lines_atom


def from_lines(lines_cell, lines_atom):
    #---------- lattice
    unit = lines_cell[0].split()[1]
    if unit[0] == '(' and unit[-1] == ')':
        unit = unit[1:-1]
    if unit.startswith('alat'):
        scale = float(unit[unit.index('=')+1:])    # in Bohr
        scale = float(Length(scale, 'bohr').to('ang'))    # in Ang
    elif unit == 'bohr':
        scale = float(Length(1.0, 'bohr').to('ang'))    # in Ang
    elif unit == 'angstrom':
        scale = 1.0    # in Ang
    else:
        ValueError('unit "{0:s}" for CELL_PARAMETERS is not supported'.format(unit))
    lattice = [[scale * float(x) for x in line.split()] for line in lines_cell[1:4]]

    #---------- species & coordinates
    unit = lines_atom[0].split()[1]
    if unit[0] == '(' and unit[-1] == ')':
        unit = unit[1:-1]
    species = []
    coords = []
    for line in lines_atom[1:]:
        fields = line.split()
        species.append(fields[0])
        coords.append([float(x) for x in fields[1:4]])
        if unit == 'crystal':
            pass    # 'coords' are already internal coordinates
        else:
            ValueError('unit "{0:s}" for ATOMIC_POSITIONS is not supported yet'.format(unit))

    structure = Structure(lattice, species, coords)

    return structure


def out_struc(fin, fout, tolerance=0.1):
    natot = get_natot(fin)
    lines_cell = extract_cell_parameters(fout)
    if lines_cell is None:    # 'relax' --> no CELL_PARAMETERS
        lines_cell = extract_cell_parameters(fin)    # get from input
    lines_atom = extract_atomic_positions(fout, natot)
    structure = from_lines(lines_cell, lines_atom)    # pymatgen format
    structure.to(fmt='poscar', filename='out_struc.vasp')
    cif = CifWriter(structure, symprec=tolerance)
    cif.write_file('out_struc.cif')


def in_struc(fin, tolerance=0.001):
    natot = get_natot(fin)
    lines_cell = extract_cell_parameters(fin)
    lines_atom = extract_atomic_positions(fin, natot)
    structure = from_lines(lines_cell, lines_atom)    # pymatgen format
    structure.to(fmt='poscar', filename='in_struc.vasp')
    cif = CifWriter(structure, symprec=tolerance)
    cif.write_file('in_struc.cif')


if __name__ == '__main__':

    if len(sys.argv) == 2:    # one argument --> structure from input file
        in_struc(sys.argv[1])
    elif len(sys.argv) == 3:    # two arguments --> structure from output file
        out_struc(sys.argv[1], sys.argv[2])
    else:
        raise SystemExit('Usage: [python] qe2vasp_cif.py pwscf.in [pwscf.out] ')
