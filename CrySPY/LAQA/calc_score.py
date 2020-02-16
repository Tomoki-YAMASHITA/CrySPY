'''
Calculate score in LAQA
'''

import numpy as np


def calc_laqa_bias(force_step, c=1.0, thres=1.0e-6):
    if force_step is None:
        return np.nan
    F1 = np.linalg.norm(force_step[-1], axis=1).mean()
    if len(force_step) == 1:
        F2 = None
    else:
        F2 = np.linalg.norm(force_step[-2], axis=1).mean()
    dF = 1.0 if F2 is None else abs(F1-F2)
    dF = thres if dF < thres else dF
    return c*F1**2/(2.0*dF)    # search maximum
