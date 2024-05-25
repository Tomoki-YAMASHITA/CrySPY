from logging import getLogger

from ..IO import io_stat, pkl_data


logger = getLogger('cryspy')

def initialize(rin):
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
