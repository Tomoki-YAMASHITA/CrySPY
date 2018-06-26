#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pymatgen import Structure


def get_spg_info(filename, tolerance=0.1):
    if filename[-5:] == '.vasp':
        with open(filename, 'r') as f:
            struc_str = f.read()
        struc = Structure.from_str(struc_str, fmt='poscar')
    else:
        struc = Structure.from_file(filename)
    spg_sym, spg_num = struc.get_space_group_info(symprec=tolerance)
    return spg_sym, spg_num


if __name__ == '__main__':
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tolerance', help='tolerance', type=float, default=0.1)
    parser.add_argument('infile', help='input file')
    args = parser.parse_args()

    # ---------- main
    filename = args.infile.split('/')[-1]    # ./aaa/bbb/POSCAR --> POSCAR
    print get_spg_info(filename, args.tolerance)
