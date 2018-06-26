#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itertools

import numpy as np
from pymatgen import Structure

from ...IO import read_input as rin


def from_file(name):
    # ---------- last structure
    with open(name, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if 'ITEM: TIMESTEP' in line:
            tmp_lines = lines[i:i+rin.natot+9]
    lines = tmp_lines

    # ---------- lattice
    xlo_bound, xhi_bound, xy = [float(x) for x in lines[5].split()]
    ylo_bound, yhi_bound, xz = [float(x) for x in lines[6].split()]
    zlo_bound, zhi_bound, yz = [float(x) for x in lines[7].split()]

    xlo = xlo_bound - min(0.0, xy, xz, xy+xz)
    xhi = xhi_bound - max(0.0, xy, xz, xy+xz)
    ylo = ylo_bound - min(0.0, yz)
    yhi = yhi_bound - max(0.0, yz)
    zlo = zlo_bound
    zhi = zhi_bound

    a = [xhi-xlo, 0, 0]
    b = [xy, yhi-ylo, 0]
    c = [xz, yz, zhi-zlo]

    lattice = np.array([a, b, c])    # in Ang, each row is a lattice vector

    # ---------- non replica coordinates
    coords = [[float(x) for x in line.split()][2:] for line in lines[9:]]

    # ---------- species
    species = [itertools.repeat(typ, times=num) for typ, num in zip(rin.atype, rin.nat)]
    species = list(itertools.chain.from_iterable(species))

    structure = Structure(lattice, species, coords)

    return structure


def write(structure, output, title="LAMMPS"):
    # ---------- convert lattice constants and angles to LAMMPS format
    cos = lambda x: np.cos(np.deg2rad(x))
    sqrt = lambda x: np.sqrt(x)

    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    lx = a
    xy = b*cos(gamma)
    xz = c*cos(beta)
    ly = sqrt(b**2 - xy**2)
    yz = (b*c*cos(alpha) - xy*xz) / ly
    lz = sqrt(c**2 - xz**2 - yz**2)

    # ---------- coordinates
    index = range(1, rin.natot+1)
    species = [itertools.repeat(typ, times=num) for typ, num in zip(index, rin.nat)]
    species = list(itertools.chain.from_iterable(species))
    sites = structure.sites

    # ---------- write
    with open(output, 'w') as f:
        f.write('# data_{0:s}\n'.format(''.join(title.split())))
        f.write('\n')
        f.write(str(rin.natot) + ' atoms\n')
        f.write(str(len(rin.atype)) + ' atom types\n')
        f.write('\n')
        f.write('0.000000  {0:10.6f}   xlo xhi\n'.format(lx))
        f.write('0.000000  {0:10.6f}   ylo yhi\n'.format(ly))
        f.write('0.000000  {0:10.6f}   zlo zhi\n'.format(lz))
        f.write('\n')
        f.write('{0:10.6f}  {1:10.6f}  {2:10.6f}   xy xz yz\n'.format(xy, xz, yz))
        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')

        for i, k, site in zip(index, species, sites):
            f.write('{0:4d}  {1:<4d}   {2[0]:7f} {2[1]:7f} {2[2]:7f}\n'.format(i, k, site.coords))
