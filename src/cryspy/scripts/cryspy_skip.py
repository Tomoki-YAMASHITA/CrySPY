#!/usr/bin/env python3
'''
Force skip if previously collected results are incorrect.
'''
import argparse
from logging import getLogger
import os

import numpy as np

from cryspy.IO import pkl_data
from cryspy.IO.out_results import out_rslt
from cryspy.start import cryspy_restart
from cryspy.util.utility import set_logger

# ---------- import later
#from cryspy.IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
#from cryspy.IO.out_results import out_laqa_energy, out_laqa_bias


def main():
    # ---------- argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('cid_list', nargs='+', type=int, help='structure ID list')
    args = parser.parse_args()

    # ---------- logger
    set_logger(
        noprint=False,
        debug=False,
        logfile='log_cryspy',
        errfile='err_cryspy',
        debugfile='debug_cryspy',
    )
    logger = getLogger('cryspy')
    logger.info('# ---------- cryspy-skip command')

    # ---------- lock
    if os.path.isfile('lock_cryspy'):
        logger.error('lock_cryspy file exists')
        raise SystemExit(1)
    else:
        with open('lock_cryspy', 'w'):
            pass    # create vacant file

    # ---------- restart
    if os.path.isfile('cryspy.stat'):
        rin, init_struc_data = cryspy_restart.restart()
    else:
        logger.error('cryspy.stat file does not exist.')
        os.remove('lock_cryspy')
        raise SystemExit(1)

    # ---------- load data and import
    logger.info('# ------ Load data')
    opt_struc_data = pkl_data.load_opt_struc()
    rslt_data = pkl_data.load_rslt()
    if rin.algo == 'RS':
        pass
    elif rin.algo == 'BO':
        opt_dscrpt_data = pkl_data.load_opt_dscrpt_data()
    elif rin.algo == 'LAQA':
        from cryspy.IO.out_results import out_laqa_status, out_laqa_step, out_laqa_score
        from cryspy.IO.out_results import out_laqa_energy, out_laqa_bias
        laqa_step = pkl_data.load_laqa_step()
        laqa_struc = pkl_data.load_laqa_struc()
        laqa_energy = pkl_data.load_laqa_energy()
        laqa_bias = pkl_data.load_laqa_bias()
        laqa_score = pkl_data.load_laqa_score()
    elif rin.algo in ['EA', 'EA-vc']:
        elite_struc = pkl_data.load_elite_struc()
        elite_fitness = pkl_data.load_elite_fitness()

    # ---------- ID loop for skip
    for cid in args.cid_list:
        logger.info(f'# ------ Skip structure ID {cid}')
        if cid not in rslt_data.index:
            logger.error(f'Structure ID {cid} does not exist in rslt_data')
            os.remove('lock_cryspy')
            raise SystemExit(1)
        # ------ common parts
        opt_struc_data[cid] = None
        rslt_data.at[cid, 'Spg_num_opt'] = 0
        rslt_data.at[cid, 'Spg_sym_opt'] = None
        rslt_data.at[cid, 'E_eV_atom'] = np.nan
        rslt_data.at[cid, 'Magmom'] = np.nan
        rslt_data.at[cid, 'Opt'] = 'skip'
        # ------ RS
        if rin.algo == 'RS':
            pass
        # ------ BO
        elif rin.algo == 'BO':
            opt_dscrpt_data[cid] = None
            pkl_data.save_opt_dscrpt_data(opt_dscrpt_data)
        # ------ LAQA
        elif rin.algo == 'LAQA':
            laqa_step[cid].append(0)
            laqa_struc[cid].append(None)
            laqa_energy[cid].append(np.nan)
            laqa_bias[cid].append(np.nan)
            laqa_score[cid].append(-float('inf'))
            pkl_data.save_laqa_step(laqa_step)
            pkl_data.save_laqa_struc(laqa_struc)
            pkl_data.save_laqa_energy(laqa_energy)
            pkl_data.save_laqa_bias(laqa_bias)
            pkl_data.save_laqa_score(laqa_score)
            out_laqa_status(laqa_step, laqa_score,
                            laqa_energy, laqa_bias)
            out_laqa_step(laqa_step)
            out_laqa_score(laqa_score)
            out_laqa_energy(laqa_energy)
            out_laqa_bias(laqa_bias)
        # ------ EA
        elif rin.algo == 'EA':
            pass
        # ------ EA-vc
        elif rin.algo == 'EA-vc':
            rslt_data.at[cid, 'Ef_eV_atom'] = np.nan
        # ------ elite
        if rin.algo in ['EA', 'EA-vc']:
            if cid in elite_struc:
                del elite_struc[cid]
                del elite_fitness[cid]
                pkl_data.save_elite_struc(elite_struc)
                pkl_data.save_elite_fitness(elite_fitness)
        # ------ common parts to save
        pkl_data.save_opt_struc(opt_struc_data)
        pkl_data.save_rslt(rslt_data)
        out_rslt(rslt_data)
        logger.info(f'Structure ID {cid} skipped successfully')

    # ---------- unlock
    os.remove('lock_cryspy')


if __name__ == '__main__':
    # ---------- main
    main()