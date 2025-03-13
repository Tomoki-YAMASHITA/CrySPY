'''
Structure file for OpenMX
written by H. Sawahata 2020/03/09
info at hikaruri.jp
'''

from pymatgen.core import Structure
from pymatgen.core.units import Length


import re


def extract_cell_parameters_from_outfile(filename):
    # ---------- last Cell vectors a1, a2, a3
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_cell = None
    for i,line in enumerate(lines):
        if re.search('a1', line):
            lines_cell = [line[7:37],lines[i+1][7:37],
                lines[i+2][7:37]]
            break
    return lines_cell


def extract_cell_parameters_from_infile(filename):
    # ---------- last Cell vectors a1, a2, a3
    # case of not optimized Cell
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_cell = None
    for i,line in enumerate(lines):
        if re.search('<Atoms.UnitVectors', line):
            lines_cell = [lines[i+1],lines[i+2],lines[i+3]]
            break
    return lines_cell

def extract_atomic_positions_from_outfile(filename, nat):
    # ---------- natot
    natot = sum(nat)
    # ---------- last atomic positions
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_atom = None
    for i, line in enumerate(lines):
        if re.search('final', line):
            ibegin = i + 4
            iend = ibegin + natot
            lines_atom = lines[ibegin:iend]
            break
    return lines_atom

def extract_atomic_positions_from_infile(filename, nat):
    # ---------- natot
    natot = sum(nat)
    # ---------- last atomic positions
    with open(filename, 'r') as f:
        lines = f.readlines()
    lines_atom = []
    for i,line in enumerate(lines):
        if re.search('<Atoms.SpeciesAndCoordinates', line):
            for j in range(0, natot):
                Atom_Row = lines[i+j+1].split()
                Row = ' '.join(Atom_Row[0:5])
                lines_atom.append(Row)
            break
    return lines_atom

def from_lines(lines_cell, lines_atom):
    # ---------- lattice
    lattice = [[float(x) for x in line.split()]
        for line in lines_cell[0:3]]
    # ---------- species & coordinates
    species = []
    coords = []
    for line in lines_atom:
        fields = line.split()
        species.append(fields[1])
        coords.append([float(x) for x in fields[2:5]])
    structure = Structure(lattice, species, coords)
    return structure

def write(rin, structure, output, mode='w'):
    # ---------- get in POSCAR format
    poscar = structure.to(fmt='poscar')
    lines = poscar.split('\n')[:-1]
    lines = [line+'\n' for line in lines]
    # ---------- write in OMX format (yet)
    with open(output, mode) as f:
        f.write('<Atoms.UnitVectors\n')
        for line in lines[2:5]:
            f.write(line)
        f.write('Atoms.UnitVectors>\n')
        f.write('<Atoms.SpeciesAndCoordinates\n')  
        for i,line in enumerate(lines[8:]):
            fields = line.split()
            up     = str(rin.upSpin[fields[3]])
            down   = str(rin.downSpin[fields[3]])
            f.write(str(i+1)+' '+fields[3]+'  '+' '.join(fields[0:3])+' '+up+' '+down+'\n')
        f.write('Atoms.SpeciesAndCoordinates>\n')  
