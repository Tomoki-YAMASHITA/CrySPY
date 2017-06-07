#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import pickle
import sys

from pymatgen import Structure
from pymatgen.io.vasp.sets import MITRelaxSet


def get_struc(filename):
    structure = Structure.from_file(filename)
    return structure


def load_init_struc(filepath):
    with open(filepath, 'rb') as dstruc:
        init_struc_data = pickle.load(dstruc)
    return init_struc_data


def kpt_check(structure, kppvol):
    mitparamset = MITRelaxSet(structure)
    kpoints = mitparamset.kpoints.automatic_density_by_vol(structure, kppvol)
    print('a =', structure.lattice.a)
    print('b =', structure.lattice.b)
    print('c =', structure.lattice.c)
    print('    Lattice vector')
    print(structure.lattice)
    print()
    print('kppvol: ', kppvol)
    print('k-points: ', kpoints.kpts[0])


if __name__ == '__main__':

    # sys.argv[1] <-- POSCAR, CONTCAR, or init_struc_data.pkl
    #                 ./aaa/bbb/POSCAR is OK
    # sys.argv[2] <-- kppvol

    nstruc = 5

    if len(sys.argv) == 3:    # two arguments
        vaspfiles = ['POSCAR', 'CONTCAR']
        filename = sys.argv[1].split('/')[-1]
        if filename in vaspfiles:
            structure = get_struc(sys.argv[1])
            kpt_check(structure, int(sys.argv[2]))
        elif filename == 'init_struc_data.pkl':
            init_struc_data = load_init_struc(sys.argv[1])
            count = 0
            for structure in init_struc_data:
                print('\n\n#---------- {}th structure'.format(count))
                kpt_check(structure, int(sys.argv[2]))
                count += 1
                if count > nstruc - 1 :
                    sys.exit()
    else:
        raise SystemExit('Usage: [python] kpt_check.py filename kppvol')
