'''
Collect results in ext
'''
from logging import getLogger
import os

from ...IO import pkl_data


logger = getLogger('cryspy')

def collect_ext():
    # ---------- get opt_struc_data and energy_data
    if os.path.isdir('ext/calc_data'):
        ext_opt_struc_data = pkl_data.load_ext_opt_struc()
        ext_energy_data = pkl_data.load_ext_energy()
    else:
        logger.error('no ext/calc_data directory')
        raise SystemExit(1)
    # ---------- return
    return ext_opt_struc_data, ext_energy_data
