'''
Restart Bayesian optimization
'''

from .select_descriptor import select_descriptor
from ..IO import pkl_data


def restart(rin, init_struc_data, prev_nstruc):
    # ---------- load BO data
    init_dscrpt_data = pkl_data.load_init_dscrpt_data()

    # ---------- get additional struc data
    tmp_struc_data = {cid: init_struc_data[cid] for cid in range(
        prev_nstruc, len(init_struc_data))}

    # ---------- calc descriptor
    tmp_dscrpt = select_descriptor(rin, tmp_struc_data)

    # ---------- update
    init_dscrpt_data.update(tmp_dscrpt)

    # ---------- save
    pkl_data.save_init_dscrpt_data(init_dscrpt_data)
