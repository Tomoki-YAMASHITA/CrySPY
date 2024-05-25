from logging import getLogger
import os

from ..IO import io_stat, pkl_data
from ..IO.out_results import out_laqa_id_hist


logger = getLogger('cryspy')


def next_selection(rin, laqa_id_data, laqa_data):
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
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)
    pkl_data.save_id_select_hist(id_select_hist)
    pkl_data.save_tot_step_select(tot_step_select)
    pkl_data.save_laqa_step(laqa_step)
    pkl_data.save_laqa_struc(laqa_struc)
    pkl_data.save_laqa_energy(laqa_energy)
    pkl_data.save_laqa_bias(laqa_bias)
    pkl_data.save_laqa_score(laqa_score)

    # ---------- status
    stat = io_stat.stat_read()
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
