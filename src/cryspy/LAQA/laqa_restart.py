'''
Restart in LAQA
'''

from logging import getLogger

from ..IO import io_stat, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def restart(stat, prev_nstruc):
    # ---------- load laqa data
    id_queueing, id_running, id_select_hist = pkl_data.load_laqa_id()
    tot_step_select, laqa_step, laqa_struc, \
        laqa_energy, laqa_bias, laqa_score = pkl_data.load_laqa_data()

    # ---------- initialize for appended ID
    for i in range(prev_nstruc, rin.tot_struc):
        laqa_step[i] = []
        laqa_struc[i] = []
        laqa_energy[i] = []
        laqa_bias[i] = []
        laqa_score[i] = [float('inf')]
        id_queueing.append(i)
    logger.info('Append scores and id_queueing')

    # ---------- status
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save for LAQA
    laqa_id_data = (id_queueing, id_running, id_select_hist)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc,
                 laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)
