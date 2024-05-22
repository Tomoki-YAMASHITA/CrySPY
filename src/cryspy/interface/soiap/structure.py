'''
Structure files in soiap
'''
import itertools

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from ...util import constants


def from_file(rin, lines, nat):
    # ---------- lattice
    lattice = [[float(x) for x in line.split()] for line in lines[1:4]]
    #     in Bohr, each column is a lattice vector
    lattice = np.array(lattice) * constants.BOHR2ANG  # Bohr --> Ang
    lattice = lattice.T    # each row is a lattice vector

    # ---------- internal coordinates
    coords = [[float(x) for x in line.split()] for line in lines[5:]]

    # ---------- species
    species = [itertools.repeat(typ, times=num) for typ, num in zip(
        rin.atype, nat)]
    species = list(itertools.chain.from_iterable(species))

    structure = Structure(lattice, species, coords)

    return structure


def write(structure, output, symprec=0.01, title="soiap"):
    # ---------- symmetrized structure
    analyzer = SpacegroupAnalyzer(structure, symprec=symprec)
    structure_sym = analyzer.get_symmetrized_structure()

    # ---------- get written data
    lattice = structure_sym.lattice
    symmetry_operations = structure_sym.spacegroup
    site_indices_inequiv = [equiv[0]
                            for equiv in structure_sym.equivalent_indices]
    sites_inequiv = [structure_sym.sites[i] for i in site_indices_inequiv]

    # ---------- write
    with open(output, 'w') as f:
        f.write('data_{0:s}\n'.format(''.join(title.split())))
        f.write('\n')
        for s, value in zip(['a', 'b', 'c'], lattice.abc):
            f.write('{0:17s} {1:12.8f}\n'.format('_cell_length_'+s, value))
        for s, value in zip(['alpha', 'beta', 'gamma'], lattice.angles):
            f.write('{0:17s} {1:12.8f}\n'.format('_cell_angle_'+s, value))
        f.write('\n')
        f.write('loop_\n')
        f.write('  _symmetry_equiv_pos_as_xyz\n')
        for symop in symmetry_operations:
            f.write("  '{0:s}'\n".format(symop.as_xyz_str()))
        f.write('\n')
        f.write('loop_\n')
        f.write('  _atom_site_type_symbol\n')
        for s in ['x', 'y', 'z']:
            f.write('  _atom_site_fract_'+s+'\n')
        for site in sites_inequiv:
            f.write('  {0:5s}  '.format(site.species_string)+' '.join(
                ['{0:9.6f}'.format(x) for x in site.frac_coords])+'\n')
