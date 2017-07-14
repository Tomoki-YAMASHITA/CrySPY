#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import math


def rndgen_lattice(spgnum, minlen, maxlen, dangle):
    #---------- for spgnum = 0: no space group
    if spgnum == 0:
        #----- trigonal = hexagonal
        crystal_system = ['Triclinic',
                          'Monoclinic',
                          'Orthorhombic',
                          'Tetragonal',
                          'Rhombohedral',
                          'Hexagonal',
                          'Cubic']
        spg = 0
        csys = random.choice(crystal_system)
    #---------- for spgnum 1--230
    else:
        spg = get_spg(spgnum)
        if 1 <= spg <= 2:
            csys = 'Triclinic'
        elif 3 <= spg <= 15:
            csys = 'Monoclinic'
        elif 16 <= spg <= 74:
            csys = 'Orthorhombic'
        elif 75 <= spg <= 142:
            csys = 'Tetragonal'
        elif 143 <= spg <= 167:
            #----- trigonal includes rhombohedral in find_wy
            csys = 'Trigonal'
        elif 168 <= spg <= 194:
            csys = 'Hexagonal'
        elif 195 <= spg <= 230:
            csys = 'Cubic'
        else:
            raise ValueError('spg is wrong')

    #---------- generate lattice constants a, b, c, alpha, beta, gamma
    if csys == 'Triclinic':
        t1 = random.uniform(minlen, maxlen)
        t2 = random.uniform(minlen, maxlen)
        t3 = random.uniform(minlen, maxlen)
        t = [t1, t2, t3]
        t.sort()
        a, b, c = t
        r = random.random()
        if r < 0.5:    # Type I
            alpha = 90.0 - random.uniform(0, dangle)
            beta  = 90.0 - random.uniform(0, dangle)
            gamma = 90.0 - random.uniform(0, dangle)
        else:    # Type II
            alpha = 90.0 + random.uniform(0, dangle)
            beta  = 90.0 + random.uniform(0, dangle)
            gamma = 90.0 + random.uniform(0, dangle)
    elif csys == 'Monoclinic':
        a = random.uniform(minlen, maxlen)
        b = random.uniform(minlen, maxlen)
        c = random.uniform(minlen, maxlen)
        if a > c:
            a, c = c, a
        alpha = gamma = 90.0
        beta = 90.0 + random.uniform(0, dangle)
    elif csys == 'Orthorhombic':
        t1 = random.uniform(minlen, maxlen)
        t2 = random.uniform(minlen, maxlen)
        t3 = random.uniform(minlen, maxlen)
        t = [t1, t2, t3]
        t.sort()
        a, b, c = t
        alpha = beta = gamma = 90.0
    elif csys == 'Tetragonal':
        a = b = random.uniform(minlen, maxlen)
        c = random.uniform(minlen, maxlen)
        alpha = beta = gamma = 90.0
    elif csys == 'Trigonal':
        a = b = random.uniform(minlen, maxlen)
        c = random.uniform(minlen, maxlen)
        alpha = beta = 90.0
        gamma = 120.0
    elif csys == 'Rhombohedral':
        a = b = c = random.uniform(minlen, maxlen)
        tangl = 90 + random.uniform(-dangle, dangle)
        alpha = beta = gamma = tangl
    elif csys == 'Hexagonal':
        a = b = random.uniform(minlen, maxlen)
        c = random.uniform(minlen, maxlen)
        alpha = beta = 90.0
        gamma = 120.0
    elif csys == 'Cubic':
        a = b = c = random.uniform(minlen, maxlen)
        alpha = beta = gamma = 90.0

    return spg, a, b, c, alpha, beta, gamma


def get_spg(spgnum):
    if spgnum == 'all':
        spg = random.randint(1, 230)
    else:
        spg = random.choice(spgnum)
    return spg


def calc_latvec(a, b, c, alpha, beta, gamma):
    #---------- degree to radian
    alpha_rad = math.radians(alpha)
    beta_rad = math.radians(beta)
    gamma_rad = math.radians(gamma)

    #---------- calculate components
    bx = b*math.cos(gamma_rad)
    by = b*math.sin(gamma_rad)
    cx = c*math.cos(beta_rad)
    cy = (c*math.cos(alpha_rad) - cx*math.cos(gamma_rad))/math.sin(gamma_rad)
    cz = math.sqrt(c*c - cx*cx - cy*cy)

    #---------- lattice vector as list
    va = [a, 0.0, 0.0]
    vb = [bx, by, 0.0]
    vc = [cx, cy, cz]

    return va, vb, vc


def calc_cos(alpha, beta, gamma):
    #---------- degree to radian
    a_rad = math.radians(alpha)
    b_rad = math.radians(beta)
    g_rad = math.radians(gamma)

    return math.cos(a_rad), math.cos(b_rad), math.cos(g_rad)
