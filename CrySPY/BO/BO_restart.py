#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import select_descriptor
from ..IO import pkl_data


def restart(init_struc_data, BO_id_data, BO_data):
    gen, next_BO_id, non_error_id, id_to_calc, id_done = BO_id_data
    descriptors, targets = BO_data

    # ---------- append non_error_id and descriptors
    non_error_id, descriptors = select_descriptor.append_X(init_struc_data, next_BO_id,
                                                           non_error_id, descriptors)
    next_BO_id = len(init_struc_data)

    # ---------- save
    BO_id_data = (gen, next_BO_id, non_error_id, id_to_calc, id_done)
    pkl_data.save_BO_id(BO_id_data)
    BO_data = (descriptors, targets)
    pkl_data.save_BO_data(BO_data)

    return BO_id_data, BO_data
