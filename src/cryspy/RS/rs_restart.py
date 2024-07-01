from ..IO import io_stat, pkl_data


def restart(rin, prev_nstruc):
    # ---------- load rs data
    id_queueing = pkl_data.load_id_queueing()

    # ---------- append id_queueing
    id_queueing.extend([i for i in range(prev_nstruc, rin.tot_struc)])

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save
    pkl_data.save_id_queueing(id_queueing)
