'''
Initialize LAQA
'''

from logging import getLogger

from ..IO import io_stat, pkl_data
from ..IO import read_input as rin


logger = getLogger('cryspy')

def initialize(stat):
    logger.info('# ---------- Initialize LAQA')

    # ---------- initialize
    tot_step_select = [0]
    laqa_step = {}
    laqa_struc = {}
    laqa_energy = {}
    laqa_bias = {}
    laqa_score = {}
    for i in range(rin.tot_struc):
        laqa_step[i] = []
        laqa_struc[i] = []
        laqa_energy[i] = []
        laqa_bias[i] = []
        laqa_score[i] = [float('inf')]
    id_queueing = [i for i in range(rin.tot_struc)]
    id_select_hist = []
    id_running = []

    # ---------- save for LAQA
    laqa_id_data = (id_queueing, id_running, id_select_hist)
    pkl_data.save_laqa_id(laqa_id_data)
    laqa_data = (tot_step_select, laqa_step, laqa_struc,
                 laqa_energy, laqa_bias, laqa_score)
    pkl_data.save_laqa_data(laqa_data)

    # ---------- status
    io_stat.set_common(stat, 'selection', 0)
    io_stat.set_common(stat, 'total_step', 0)
    io_stat.set_id(stat, 'selected_id', id_queueing)    # all IDs
    io_stat.set_id(stat, 'id_queueing', id_queueing)    # all IDs
    io_stat.write_stat(stat)

    # ---------- log
    logger.info('# ---------- Selection 0')
    if len(id_queueing) > 30:
        logger.info(f'selected_id: {len(id_queueing)} IDs')
    else:
        x = ' '.join(str(a) for a in id_queueing)
        logger.info(f'selected_id: {x}')
