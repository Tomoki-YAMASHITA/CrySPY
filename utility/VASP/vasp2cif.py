#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from pymatgen import Structure
from pymatgen.io.cif import CifWriter


def get_cif(filename, tolerance=0.1):
    opt_struc = Structure.from_file(filename)
    cif = CifWriter(opt_struc, symprec=tolerance)
    cif.write_file('struc.cif')


if __name__ == '__main__':

    if len(sys.argv) == 2:    # only one argument
        get_cif(sys.argv[1])
    elif len(sys.argv) == 3:    # two arguments --> change tolerance
        tolerance = float(sys.argv[2])
        get_cif(sys.argv[1], tolerance)
    else:
        raise SystemExit('Usage: [python] vasp2cif.py filename [tolerance]')
