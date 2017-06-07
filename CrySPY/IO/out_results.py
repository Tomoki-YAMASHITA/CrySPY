#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import read_input as rin


def write_rslt(rslt_data):
    #---------- asc in Struc_ID or Gen
    with open('./data/cryspy_rslt', 'w') as frslt:
        if rin.algo == 'RS':
            frslt.write(rslt_data.sort_values(by=['Struc_ID'], ascending=True).to_string(index=False))
        elif rin.algo == 'BO':
            frslt.write(rslt_data.sort_values(by=['Gen', 'Struc_ID'], ascending=True).to_string(index=False))

    #---------- asc in energy
    with open('./data/cryspy_rslt_energy_asc', 'w') as fasc:
        fasc.write(rslt_data.sort_values(by=['Energy'], ascending=True).to_string(index=False))
