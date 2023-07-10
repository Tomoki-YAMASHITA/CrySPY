'''
Selection in LAQA
'''

from logging import getLogger
import os

from ..IO import io_stat, pkl_data
from ..IO import read_input as rin
from ..IO.out_results import out_laqa_id_hist


logger = getLogger('cryspy')

def next_selection(stat, laqa_id_data, laqa_data):
    # ---------- laqa_id_data and laqa_data
    (id_queueing, id_running, id_select_hist) = laqa_id_data
    (tot_step_select, laqa_step, laqa_struc,
     laqa_energy, laqa_bias, laqa_score) = laqa_data

    # ---------- LAQA selection
    for k, v in sorted(laqa_score.items(), key=lambda x: -x[1][-1]):
        if v[-1] == -float('inf'):
            break
        else:
            id_queueing.append(k)
            if len(id_queueing) == rin.nselect_laqa:
                break

    # ---------- done LAQA
    if len(id_queueing) == 0:
        logger.info('\nDone all structures!')
        os.remove('lock_cryspy')
        raise SystemExit()

    # ---------- append id_select_hist and out
    id_select_hist.append(id_queueing[:])    # append shallow copy
    out_laqa_id_hist(id_select_hist)

    # ---------- tot_step_select for next selection
    tot_step_select.append(0)

    # ---------- save
    laqa_id_data = (id_queueing, id_running, id_select_hist)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc,
                 laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)

    # ---------- status
    io_stat.set_common(stat, 'selection', len(id_select_hist))
    io_stat.set_id(stat, 'selected_id', id_queueing)
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- out and log
    logger.info(f'# ---------- Selection {len(id_select_hist)}')
    if len(id_queueing) > 30:
        logger.info(f'selected_id: {len(id_queueing)} IDs')
    else:
        x = ' '.join(str(a) for a in id_queueing)
        logger.info(f'selected_id: {x}')
