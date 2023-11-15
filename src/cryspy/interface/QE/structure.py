'''
Structure file for Quantum ESPRESSO
'''

from logging import getLogger

from pymatgen.core import Structure
from pymatgen.core.units import Length

from ...IO import read_input as rin


logger = getLogger('cryspy')

def extract_cell_parameters(filename):
    # ---------- last CELL_PARAMETERS
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


def extract_atomic_positions(filename, nat):
    # ---------- last ATOMIC_POSITIONS
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_atom = None
    natot = sum(nat)    # do not use rin.natot here for EA-vc
    for i, line in enumerate(reversed(lines)):
        if 'ATOMIC_POSITIONS' in line:
            ibegin = len(lines) - 1 - i
            iend = ibegin + natot + 1
            lines_atom = lines[ibegin:iend]
            break
    return lines_atom


def from_lines(lines_cell, lines_atom):
    # ---------- lattice
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
        logger.error(f'unit "{unit:s}" for CELL_PARAMETERS'
                         ' is not supported')
    lattice = [[scale * float(x) for x in line.split()]
               for line in lines_cell[1:4]]

    # ---------- species & coordinates
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
            logger.error(f'unit "{unit:s}" for ATOMIC_POSITIONS'
                             ' is not supported yet')

    structure = Structure(lattice, species, coords)
    return structure


def write(structure, output, mode='w'):
    # ---------- get in POSCAR format
    poscar = structure.to(fmt='poscar')
    lines = poscar.split('\n')[:-1]
    lines = [line+'\n' for line in lines]

    # ---------- write in QE format
    with open(output, mode) as f:
        f.write('CELL_PARAMETERS angstrom\n')
        for line in lines[2:5]:
            f.write(line)
        f.write('ATOMIC_POSITIONS crystal\n')
        for line in lines[8:]:
            fields = line.split()
            f.write(fields[3]+'  '+' '.join(fields[0:3])+'\n')
