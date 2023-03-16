'''
Restart Bayesian optimization
'''

from .select_descriptor import select_descriptor
from ..IO import pkl_data


def restart(init_struc_data, prev_nstruc):
    # ---------- load BO data
    (init_dscrpt_data, opt_dscrpt_data,
     bo_mean, bo_var, bo_score) = pkl_data.load_bo_data()

    # ---------- get additional struc data
    tmp_struc_data = {cid: init_struc_data[cid] for cid in range(
        prev_nstruc, len(init_struc_data))}

    # ---------- calc descriptor
    tmp_dscrpt = select_descriptor(tmp_struc_data)

    # ---------- update
    init_dscrpt_data.update(tmp_dscrpt)

    # ---------- save
    bo_data = (init_dscrpt_data, opt_dscrpt_data, bo_mean, bo_var, bo_score)
    pkl_data.save_bo_data(bo_data)
