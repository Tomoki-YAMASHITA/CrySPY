from ..IO import io_stat, pkl_data


def restart(rin, prev_nstruc):
    # ---------- load rs data
    id_queueing, id_running = pkl_data.load_rs_id()

    # ---------- append id_queueing
    id_queueing.extend([i for i in range(prev_nstruc, rin.tot_struc)])

    # ---------- status
    stat = io_stat.stat_read()
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save
    rs_id_data = (id_queueing, id_running)
    pkl_data.save_rs_id(rs_id_data)

    # ---------- ext
    if rin.calc_code == 'ext':
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('out\n')
