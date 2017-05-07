#!/usr/bin/env python
# -*- coding: utf-8 -*-


def write_poscar(va, vb, vc, incoord, atype, nat):
    #---------- float --> character
    va_c = '    {: .10f}      {: .10f}      {: .10f}\n'.format(va[0], va[1], va[2])
    vb_c = '    {: .10f}      {: .10f}      {: .10f}\n'.format(vb[0], vb[1], vb[2])
    vc_c = '    {: .10f}      {: .10f}      {: .10f}\n'.format(vc[0], vc[1], vc[2])

    #---------- write
    with open('POSCAR', 'w') as poscar:
        poscar.write('cpsy\n')
        poscar.write('1.0\n')
        poscar.write(va_c)
        poscar.write(vb_c)
        poscar.write(vc_c)
        for i in atype:
            poscar.write('{:>8}'.format(i))
        poscar.write('\n')
        for i in nat:
            poscar.write('{:>8d}'.format(i))
        poscar.write('\n')
        poscar.write('Direct\n')
        for j in incoord:    # atom loop
            for i in j:
                poscar.write('    {: .10f}'.format(i))
            poscar.write('\n')


