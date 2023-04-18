'''
Calculate score in LAQA
'''

import numpy as np


def calc_laqa_bias(force_step, stress_step, wf=0.1, ws=10.0, thres=1.0e-6):
    '''
    force_step: force_step_data at a certain ID, selection (stage)
    stress_step: stress_step_data at a certain ID, selection (stage)

    In this function,
    force_step[step][atom]
    stress_step[step][atom]
    '''
    # ---------- force
    if force_step is None:
        return np.nan
    F1 = np.linalg.norm(force_step[-1], axis=1).mean()
    if len(force_step) == 1:
        F2 = None
    else:
        F2 = np.linalg.norm(force_step[-2], axis=1).mean()
    dF = 1.0 if F2 is None else abs(F1-F2)
    dF = thres if dF < thres else dF

    # ---------- stress
    if stress_step is None:
        return np.nan
    S = np.abs(stress_step[-1]).mean()    # use abs. valule

    # ------ return score
    return wf*F1**2/(2.0*dF) + ws*S    # search maximum