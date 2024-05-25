from logging import getLogger

from ..IO import io_stat, pkl_data


logger = getLogger('cryspy')

def restart(rin, prev_nstruc):
    # ---------- load laqa data
    id_queueing = pkl_data.load_id_queueing()
    laqa_step = pkl_data.load_laqa_step()
    laqa_struc = pkl_data.load_laqa_struc()
    laqa_energy = pkl_data.load_laqa_energy()
    laqa_bias = pkl_data.load_laqa_bias()
    laqa_score = pkl_data.load_laqa_score()

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
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save for LAQA
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_laqa_step(laqa_step)
    pkl_data.save_laqa_struc(laqa_struc)
    pkl_data.save_laqa_energy(laqa_energy)
    pkl_data.save_laqa_bias(laqa_bias)
    pkl_data.save_laqa_score(laqa_score)
