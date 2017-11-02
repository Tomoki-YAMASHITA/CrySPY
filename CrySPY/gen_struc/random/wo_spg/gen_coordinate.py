#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import random


def rndgen_coord(natot, va, vb, vc, mindist, maxcnt):
    # ---------- generate internal coordinates
    cnt = 0
    incoord = []
    while len(incoord) < natot:
        tmp_coord = [random.random() for i in range(3)]
        if not incoord:
            incoord.append(tmp_coord)    # first atom automatically registered
        else:
            cnt = check_distance(va, vb, vc, incoord, tmp_coord, cnt, mindist)
            if cnt == 0:
                incoord.append(tmp_coord)
        if maxcnt < cnt:
            return []
    return incoord


def check_distance(va, vb, vc, incoord, tmp_coord, cnt, mindist):
    for j in incoord:    # atom loop
        # ---------- check interatomic distance
        dist = calc_atom_dist(va, vb, vc, j, tmp_coord)
        if dist < mindist:
            return cnt + 1
    return 0


def calc_atom_dist(va, vb, vc, incoordA, incoordB):
    # ---------- search nearest neighbor
    din = []
    for i in range(3):
        dtmp = incoordA[i] - incoordB[i]
        if dtmp < -0.5:
            dtmp += 1.0
        elif 0.5 < dtmp:
            dtmp -= 1.0
        din.append(dtmp)

    # ---------- calculate distance
    dx = din[0]*va[0] + din[1]*vb[0] + din[2]*vc[0]
    dy = din[0]*va[1] + din[1]*vb[1] + din[2]*vc[1]
    dz = din[0]*va[2] + din[1]*vb[2] + din[2]*vc[2]
    dist = math.sqrt(dx*dx + dy*dy + dz*dz)

    return dist
