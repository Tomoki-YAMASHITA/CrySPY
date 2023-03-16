'''
Restart random search
'''

from ..IO import io_stat, pkl_data
from ..IO import read_input as rin


def restart(stat, prev_nstruc):
    # ---------- load rs data
    id_queueing, id_running = pkl_data.load_rs_id()

    # ---------- append id_queueing
    id_queueing.extend([i for i in range(prev_nstruc, rin.tot_struc)])

    # ---------- status
    io_stat.set_id(stat, 'id_queueing', id_queueing)
    io_stat.write_stat(stat)

    # ---------- save
    rs_id_data = (id_queueing, id_running)
    pkl_data.save_rs_id(rs_id_data)

    # ---------- ext
    if rin.calc_code == 'ext':
        with open('ext/stat_job', 'w') as fstat:
            fstat.write('out\n')
