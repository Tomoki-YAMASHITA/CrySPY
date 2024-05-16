from ..IO import io_stat, pkl_data


def initialize(rin):
    # ---------- initialize
    id_queueing = [i for i in range(rin.tot_struc)]
    id_running = []

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save
    rs_id_data = (id_queueing, id_running)
    pkl_data.save_rs_id(rs_id_data)
