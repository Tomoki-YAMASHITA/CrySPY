'''
Collect results in ext
'''

import os

from ...IO import pkl_data


def collect_ext():
    # ---------- get opt_struc_data and energy_data
    if os.path.isdir('ext/calc_data'):
        ext_opt_struc_data = pkl_data.load_ext_opt_struc()
        ext_energy_data = pkl_data.load_ext_energy()
    else:
        raise IOError('no ext/calc_data directory')

    # ---------- return
    return ext_opt_struc_data, ext_energy_data
