#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen import Structure
from pymatgen.io.cif import CifWriter


def get_cif(filename, tolerance=0.1):
    struc = Structure.from_file(filename)
    cif = CifWriter(struc, symprec=tolerance)
    cif.write_file(filename+'.cif')


if __name__ == '__main__':
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tolerance', help='tolerance', type=float, default=0.1)
    parser.add_argument('infile', help='input file')
    args = parser.parse_args()

    # ---------- main
    filename = args.infile.split('/')[-1]    # ./aaa/bbb/POSCAR --> POSCAR
    get_cif(filename, args.tolerance)
