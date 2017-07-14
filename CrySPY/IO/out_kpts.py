#!/usr/bin/env python
# -*- coding: utf-8 -*-


def write_kpts(kpt_data):
    #---------- asc in ID
    with open('./data/kpts_rslt', 'w') as frslt:
        frslt.write('{0:>10}  {1:>10}\n'.format('Struc_ID', 'k-points'))
        for key, value in sorted(kpt_data.items()):
            frslt.write('{0:>10}  {1}\n'.format(key, value))
