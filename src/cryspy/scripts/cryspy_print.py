#!/usr/bin/env python3
import argparse
from dataclasses import fields
import gzip
from logging import getLogger
from math import gcd
from pathlib import Path
import pickle
from pprint import pprint

from cryspy.IO.out_results import out_rslt
from cryspy.util.utility import set_logger


PPRINT_FILES = {
    'id_queueing.pkl',
    'id_running.pkl',
    'energy_step_data.pkl',
    'force_step_data.pkl',
    'stress_step_data.pkl',
    'kpt_data.pkl',
    'rng_state_data.pkl',
    # ------ EA
    'gen.pkl',
    'ea_info.pkl',
    'ea_origin.pkl',
    'elite_fitness.pkl',
    # ------ EA-vc
    'nat_data.pkl',
    # ------ BO
    'n_selection.pkl',
    'id_select_hist.pkl',
    'bo_mean.pkl',
    'bo_var.pkl',
    'bo_score.pkl',
    # ------ LAQA
    'tot_step_select.pkl',
    'laqa_energy.pkl',
    'laqa_bias.pkl',
    'laqa_score.pkl',
    'laqa_step.pkl',
    'laqa_struc.pkl',
}

STRUCTURE_FILES = {
    'init_struc_data.pkl',
    'opt_struc_data.pkl',
    'struc_step_data.pkl',
    'elite_struc.pkl',
    'init_dscrpt_data.pkl',
    'opt_dscrpt_data.pkl',
}

INPUT_FILES = {
    'input_data.pkl',
}

RESULT_FILES = {
    'rslt_data.pkl',
}

HDIST_FILES = {
    'hdist_data.pkl',
}

PD_FILES = {
    'pd_data.pkl',
}

CN_COMB_FILES = {
    'cn_comb_data.pkl',
}


def extract_pkl_name(filepath):
    path = Path(filepath)
    filename = path.name
    if filename.endswith('.gz'):
        filename = Path(filename).stem
    return filename


def out_input(rin):
    print('[basic]')
    for data_field in fields(rin):
        key = data_field.name
        if key == 'struc_mode':
            print('')
            print('[structure]')
        if key == 'ymax':
            print('')
            print('[visual]')
        if key == 'check_mindist_opt':
            print('')
            print('[option]')
        if rin.algo == 'BO' and key == 'nselect_bo':
            print('')
            print('[BO]')
        if rin.algo == 'LAQA' and key == 'nselect_laqa':
            print('')
            print('[LAQA]')
        if rin.algo in ['EA', 'EA-vc'] and key == 'n_pop':
            print('')
            print('[EA]')
        if rin.calc_code == 'VASP' and key == 'kpt_flag':
            print('')
            print('[VASP]')
        if rin.calc_code == 'QE' and key == 'kpt_flag':
            print('')
            print('[QE]')
        if rin.calc_code == 'OMX' and key == 'kpt_flag':
            print('')
            print('[OMX]')
        if rin.calc_code == 'soiap' and key == 'kpt_flag':
            print('')
            print('[soiap]')
        if rin.calc_code == 'LAMMPS' and key == 'kpt_flag':
            print('')
            print('[LAMMPS]')
        if rin.calc_code == 'ASE' and key == 'kpt_flag':
            print('')
            print('[ASE]')
        if rin.calc_code == 'ext' and key == 'kpt_flag':
            print('')
            print('[ext]')
        if getattr(rin, key) is not None:
            print(f'{key} = {getattr(rin, key)}')


def _reduced_composition(composition):
    divisor = gcd(*composition)
    return tuple(value // divisor for value in composition)


def print_rslt_data(
        rslt_data,
        no_sort=False,
        head=None,
        tail=None,
        all_rows=False,
        composition=None):
    # ---------- constant
    DEFAULT_MAX_ROWS = 100

    # ---------- composition filter
    if composition is None or rslt_data.empty:
        df = rslt_data
    else:
        target = _reduced_composition(composition)
        comp_filter = rslt_data['Num_atom'].apply(
            lambda nat: _reduced_composition(nat) == target
        )
        df = rslt_data[comp_filter]

    # ---------- options
    if no_sort:
        sort_msg = 'Not sorted'
    else:
        order = 'Ef_eV_atom' if 'Ef_eV_atom' in df.columns else 'E_eV_atom'
        df = df.sort_values(by=[order], ascending=True)
        sort_msg = f'Sorted by: {order}'
    if composition is not None and df.empty:
        print(f'No results found for composition: {tuple(composition)}')
        return
    if head is not None:
        print(df.head(head).to_string())
        return
    if tail is not None:
        print(df.tail(tail).to_string())
        return
    if all_rows or len(df) <= DEFAULT_MAX_ROWS:
        print(df.to_string())
        return

    # ---------- print first DEFAULT_MAX_ROWS rows
    print(f'Number of structures: {len(df)}')
    print(sort_msg)
    print(f'Showing first {DEFAULT_MAX_ROWS} rows. Use --head N, --tail N, or --all to display more.')
    print('')
    print(df.head(DEFAULT_MAX_ROWS).to_string())


def main():
    '''
    print pkl_data
    '''
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='input file')
    parser.add_argument('--no-sort', action='store_true',
                        help='do not sort rslt_data.pkl')
    row_selection = parser.add_mutually_exclusive_group()
    row_selection.add_argument('--head', type=int,
                               help='show first N rows of rslt_data.pkl')
    row_selection.add_argument('--tail', type=int,
                               help='show last N rows of rslt_data.pkl')
    row_selection.add_argument('--all', action='store_true',
                               help='show all rows of rslt_data.pkl')
    parser.add_argument(
        '-c', '--composition',
        type=int,
        nargs='+',
        help='composition ratio in atype order for EA-vc results, e.g. -c 1 2',
    )
    parser.add_argument('--write', action='store_true',
                        help='write cryspy_rslt files from rslt_data.pkl')
    args = parser.parse_args()
    # ------ check options for rslt_data.pkl
    if args.head is not None and args.head <= 0:
        parser.error('--head must be a positive integer')
    if args.tail is not None and args.tail <= 0:
        parser.error('--tail must be a positive integer')
    if args.write and any([
            args.no_sort,
            args.head is not None,
            args.tail is not None,
            args.all,
            args.composition is not None]):
        parser.error('--write cannot be used with display options')
    if args.composition is not None:
        if any(value < 0 for value in args.composition):
            parser.error('--composition values must be non-negative integers')
        if not any(args.composition):
            parser.error('--composition must contain at least one positive integer')
    if args.write and not Path('cryspy.in').is_file():
        parser.error(
            '--write must be run in a directory containing cryspy.in'
        )

    # ---------- logger for --write
    if args.write:
        set_logger(
            noprint=False,
            debug=False,
            logfile='log_cryspy',
            errfile='err_cryspy',
            debugfile='debug_cryspy',
        )
        logger = getLogger('cryspy')
        logger.info('# ---------- cryspy-print command with --write')

    # ---------- extract pkl_name
    #     e.g.
    #     ./data/pkl_data/init_struc_data.pkl --> init_struc_data.pkl
    #     ./data/pkl_data/init_struc_data.pkl.gz --> init_struc_data.pkl
    pkl_name = extract_pkl_name(args.infile)

    # ---------- load pkl data
    if args.infile.endswith('.gz'):
        with gzip.open(args.infile, 'rb') as f:
            pkl_data = pickle.load(f)
    else:
        with open(args.infile, 'rb') as f:
            pkl_data = pickle.load(f)
    if args.composition is not None and pkl_name not in RESULT_FILES:
        parser.error('--composition can be used only with rslt_data.pkl')

    # ---------- print pkl data
    if pkl_name in PPRINT_FILES:
        pprint(pkl_data)

    elif pkl_name in STRUCTURE_FILES:
        print(f'Number of structures: {len(pkl_data)}')

    elif pkl_name in INPUT_FILES:
        out_input(pkl_data)

    elif pkl_name in RESULT_FILES:
        if args.composition is not None:
            if 'Num_atom' not in pkl_data.columns:
                parser.error('--composition is available only for EA-vc results')
            num_atom_data = pkl_data['Num_atom'].dropna()
            if (
                    not num_atom_data.empty
                    and len(args.composition) != len(num_atom_data.iloc[0])):
                parser.error(
                    '--composition must have the same number of values as atype'
                )
        if args.write:
            out_rslt(pkl_data)
            logger.info('Generated ./data/cryspy_rslt')
            logger.info('Generated ./data/cryspy_rslt_energy_asc')
        else:
            print_rslt_data(
                pkl_data,
                no_sort=args.no_sort,
                head=args.head,
                tail=args.tail,
                all_rows=args.all,
                composition=args.composition)

    elif pkl_name in HDIST_FILES:
        if pkl_data:
            latest_gen = max(pkl_data.keys())
            hdist = pkl_data[latest_gen]
            print(f'Latest generation: {latest_gen}')
            print(f'Number of structures: {len(hdist)}')
        else:
            print('No hull distance data')

    elif pkl_name in PD_FILES:
        if pkl_data:
            latest_gen = max(pkl_data.keys())
            phase_diagram = pkl_data[latest_gen]
            print(f'Latest generation: {latest_gen}')
            print(f'Number of stable entries: {len(phase_diagram.stable_entries)}')
            print(f'Number of all entries: {len(phase_diagram.all_entries)}')
        else:
            print('No phase diagram data')

    elif pkl_name in CN_COMB_FILES:
        print(f'mode = {pkl_data["mode"]}')
        print(f'cn_mode = {pkl_data["cn_mode"]}')
        print(f'll_nat = {pkl_data["ll_nat"]}')
        print(f'ul_nat = {pkl_data["ul_nat"]}')
        print(f'charge = {pkl_data["charge"]}')
        print(f'max_cn_grid_points = {pkl_data["max_cn_grid_points"]}')
        if pkl_data['cn_comb'] is None:
            print('cn_comb = None')
        else:
            print(f'Number of charge-neutral combinations: {len(pkl_data["cn_comb"])}')

    else:
        raise ValueError(f'Unsupported pkl file: {pkl_name}')


if __name__ == '__main__':
    main()
