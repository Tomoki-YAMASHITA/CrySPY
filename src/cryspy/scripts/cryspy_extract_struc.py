#!/usr/bin/env python3
#
# (old) extract_struc.py --> cryspy_extract_struc.py
#
#   2026/05/16 T. Yamashita
#   Integrated into CrySPY as the subcommand `cryspy-extract-struc`
#
#   2024/04/16 T. Yamashita
#   --tolerance option
#   gzip file
#
#   2023/07/21 T. Yamashita
#   --print option
#
#   2023/07/19 T. Yamashita
#   For pandas, pkl_load() --> pd.read_pickle()
#
#   2022/12/01 T. Yamashita
#   --top option
#
#   2022/05/15 T. Yamashita
#   --all option
#
#   2020/05/26 T. Yamashita
#   extract a structure from init_struc_data.pkl or opt_struc_data.pkl
#   pymatgen is required
#
import argparse
import gzip
import os
import pickle
import sys

import pandas as pd


def main():
    '''
    Extract structures from init_struc_data.pkl or opt_struc_data.pkl.

    example:
    - print ID 7 10 12
      cryspy-extract-struc opt_struc_data.pkl -i 7 10 12 -p

    - write cifs of ID 7 10 12 (output 7.cif, 10.cif, 12.cif)
      cryspy-extract-struc opt_struc_data.pkl -i 7 10 12
    - write cifs of ID 7 10 12 (output 7.cif, 10.cif, 12.cif) with symmetry information
      cryspy-extract-struc opt_struc_data.pkl -i 7 10 12 -s

    top k needs rslt_data.pkl
    - write cif of top k structures (k = 3)
      cryspy-extract-struc opt_struc_data.pkl -t 3
    - write cif of top k structures (k = 3) with the rank in the filename
      cryspy-extract-struc opt_struc_data.pkl -t 3 -r
    - write cif of top k structures (k = 3) with the rank in the filename and symmetry information
      cryspy-extract-struc opt_struc_data.pkl -t 3 -rs

    - write all (output 0.cif, 1.cif, 2.cif ....)
      cryspy-extract-struc opt_struc_data.pkl -a
    - write all (output 0.cif, 1.cif, 2.cif ....) with symmetry information
      cryspy-extract-struc opt_struc_data.pkl -as
    '''
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-p', '--print',
        help='just print, e.g., cryspy-extract-struc opt_struc_data.pkl -i 7 10 12 -ps',
        action='store_true',
    )
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument(
        '-a', '--all_id',
        help='all structures, e.g., cryspy-extract-struc opt_struc_data.pkl -as',
        action='store_true',
    )
    selection.add_argument(
        '-i', '--index',
        help='structure ID, e.g., cryspy-extract-struc opt_struc_data.pkl -i 7 10 12 -s',
        type=int,
        nargs='+',
    )
    selection.add_argument(
        '-t', '--top',
        help='top k structures, e.g. (k = 3), cryspy-extract-struc opt_struc_data.pkl -t 3 -s',
        type=int,
    )
    parser.add_argument(
        '-r', '--rank',
        help='add rank in file names, e.g., cryspy-extract-struc opt_struc_data.pkl -t 3 -rs',
        action='store_true',
    )
    parser.add_argument(
        '-s', '--symmetrized',
        help='symmetrized structure, e.g., cryspy-extract-struc opt_struc_data.pkl -i 7 10 12 -s',
        action='store_true',
    )
    parser.add_argument(
        '--tolerance',
        help=(
            'tolerance for symmetrization (default 0.01), '
            'e.g., cryspy-extract-struc opt_struc_data.pkl -i 0 1 -s --tolerance 0.01'
        ),
        type=float,
        default=0.01,
    )
    parser.add_argument('infile', help='input file')
    args = parser.parse_args()

    # ---------- load struc_data
    if not os.path.isfile(args.infile):
        parser.error(f'{args.infile} does not exist')
    if args.infile.endswith('.gz'):
        with gzip.open(args.infile, 'rb') as f:
            struc_data = pickle.load(f)
    else:
        with open(args.infile, 'rb') as f:
            struc_data = pickle.load(f)

    # ---------- index
    if args.index:   # not vacant
        for cid in args.index:
            if cid not in struc_data:
                print(f'Warning: ID {cid} does not exist. Skipped.', file=sys.stderr)
                continue
            if args.print:
                print(f'\nID {cid}')
                print(struc_data[cid])
            elif args.symmetrized:
                struc_data[cid].to(fmt='cif', filename=f'{cid}.cif', symprec=args.tolerance)
            else:
                struc_data[cid].to(fmt='cif', filename=f'{cid}.cif')
        raise SystemExit()

    # ---------- top k
    if args.top:   # not vacant
        # ------ load rslt_data.pkl. It must be in the same directory as the input
        rslt_path = os.path.join(os.path.dirname(os.path.abspath(args.infile)), 'rslt_data.pkl')
        if not os.path.isfile(rslt_path):
            parser.error(f'{rslt_path} does not exist')
        rslt_data = pd.read_pickle(rslt_path)
        # ------ top k data
        top_k = rslt_data.sort_values(by=['E_eV_atom'], ascending=True).head(args.top)
        for k, cid in enumerate(top_k.index):
            if args.print:
                print(f'\nID {cid}')
                print(struc_data[cid])
            else:
                if args.rank:
                    cifname = f'{k+1}_{cid}.cif'
                else:
                    cifname = f'{cid}.cif'
                if args.symmetrized:
                    struc_data[cid].to(fmt='cif', filename=cifname, symprec=args.tolerance)
                else:
                    struc_data[cid].to(fmt='cif', filename=cifname)
        raise SystemExit()

    # ---------- all
    if args.all_id:
        for cid, struc in struc_data.items():
            if args.print:
                print(f'\nID {cid}')
                print(struc_data[cid])
            elif args.symmetrized:
                struc.to(fmt='cif', filename=f'{cid}.cif', symprec=args.tolerance)
            else:
                struc.to(fmt='cif', filename=f'{cid}.cif')


if __name__ == '__main__':
    main()