from ..IO import io_stat, pkl_data


def initialize(rin):
    # ---------- initialize
    id_queueing = [i for i in range(1, rin.tot_struc + 1)]
    id_running = []

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save
    pkl_data.save_id_queueing(id_queueing)
    pkl_data.save_id_running(id_running)
