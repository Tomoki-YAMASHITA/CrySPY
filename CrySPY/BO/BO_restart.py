#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import select_descriptor
from ..IO import pkl_data


def restart(init_struc_data, BO_id_data, BO_data, prev_nstruc):
    gen, non_error_id, id_to_calc, id_done = BO_id_data
    descriptors, targets = BO_data

    # ---------- append non_error_id and descriptors
    non_error_id, descriptors = select_descriptor.append_X(init_struc_data, prev_nstruc,
                                                           non_error_id, descriptors)

    # ---------- save
    BO_id_data = (gen, non_error_id, id_to_calc, id_done)
    pkl_data.save_BO_id(BO_id_data)
    BO_data = (descriptors, targets)
    pkl_data.save_BO_data(BO_data)

    return BO_id_data, BO_data
