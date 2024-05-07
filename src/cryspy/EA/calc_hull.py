from logging import getLogger

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

from ..IO import read_input as rin


logger = getLogger('cryspy')


def calc_convex_hull_2d(ratio_data, ef_all, c_ids, gen):
    '''
    Input:
        ratio_data [dict]: ratio of all structures, {ID: [ratio list], ...}
        ef_all [dict]: formation energy of all structures, {ID: Ef, ...}
        c_ids [array]: ID array of current generation structures
        gen [int]: current generation

    Return:
        hdist [dict]: hull distance of all structures, {ID: distance, ...}
    '''
    # ---------- initialize
    ratio_chull =  [0.0, 1.0]    # only ef < 0 for calculation of convex hull
    ef_chull =  [0.0, 0.0]       # only ef < 0 for calculation of convex hull

    # ---------- dict to list
    for cid in ratio_data:
        if ef_all[cid] < 0.0:
            ratio_chull.append(ratio_data[cid][0])
            ef_chull.append(ef_all[cid])

    # ---------- no negative Ef
    #            in this case, hull distance is equivalent to Ef
    if len(ef_chull) == 2:    # only end points
        # np.nan --> np.inf in hdist
        hdist = {cid: np.inf if np.isnan(ef_all[cid]) else ef_all[cid] for cid in ef_all}
        draw_convex_hull_2d(None, ratio_data, ef_all, c_ids, gen)
        return hdist

    # ---------- calc convex hull
    points = list(zip(ratio_chull, ef_chull))
    hull = ConvexHull(points)
    vpoints = hull.points[hull.vertices]    # hull.vertices: index of vertices in hull.points
    vpoints = np.vstack((vpoints, vpoints[0]))    # just for plot. vpoints[0] should be [1, 0] by ConvexHull()

    # ---------- hull distance
    hdist = {
        cid: np.inf if np.isnan(ef_all[cid])
        else hull_distance_2d(ratio_data[cid][0], ef_all[cid], hull.equations)
        for cid in ef_all
    }

    # ---------- draw convex hull
    draw_convex_hull_2d(vpoints, ratio_data, ef_all, c_ids, gen)

    # ---------- return
    return hdist


def hull_distance_2d(x, y, equations):
    '''
    equations: [eq0, eq1, eq2, ...], eq0: [a, b, c] for a*x + b*y + c = 0

    Find the distance from all equations and adopt the minimum distance.
    '''
    hdists = []
    for eq in equations:
        if eq[0] != 0.0:
            dist = y - (-eq[0]*x - eq[2])/eq[1]
            if dist < 0.0:
                logger.warning('hdist <= 0.0, check hull distance.')
            hdists.append(dist)

    return min(hdists)


def draw_convex_hull_2d(vpoints, ratio_data, ef_all, c_ids, gen):
    # ---------- setting
    plt.rcParams.update(_set_params())

    # ---------- fig
    fig, ax = plt.subplots(1, 1)

    # ---------- hline
    ax.axhline(y=0, xmin=0, xmax=1, color='black', linestyle='--')

    # ---------- label
    ax.set_xlabel('$x$ in '+f'{rin.atype[0]}'+'$_{x}$'+f'{rin.atype[1]}'+'$_{1-x}$')
    ax.set_ylabel('Formation energy (eV/atom)')

    # ---------- lim
    min_ef = min(ef_all.values())
    ymin = min_ef * 1.1 if min_ef < -0.01 else -0.01
    ax.set_ylim(ymin, 0.05)

    # ---------- plot ef_all
    for cid in ef_all:
        if np.isnan(ef_all[cid]):
            continue
        if cid in c_ids:
            ax.plot(ratio_data[cid][0], ef_all[cid], 'o', ms=12, mew=2.0, c='C2', alpha=0.8)
        else:
            ax.plot(ratio_data[cid][0], ef_all[cid], 'o', ms=12, mew=2.0, c='C0', alpha=0.8)

    # ---------- plot convex hull
    if vpoints is None:    # no negative Ef, plot only end points
        ax.plot([0, 1], [0, 0], 'o', ms=12, mew=2.0, c='C1')
    else:
        ax.plot(vpoints[1:, 0], vpoints[1:, 1], '-o', ms=12, mew=2.0, c='C1')

    # ---------- save figure
    fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.png', bbox_inches='tight')


def _set_params():
    rcParams_dict = {
        # ---------- figure
        'figure.figsize': [8, 6],
        'figure.dpi': 120,
        'figure.facecolor': 'white',
        # ---------- axes
        'axes.grid': True,
        'axes.linewidth': 1.5,
        # ---------- ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 8.0,
        'ytick.major.size': 8.0,
        # ---------- lines
        'lines.linewidth': 2.5,
        'lines.markersize': 12,
        # ---------- grid
        'grid.linestyle': ':',
        # ---------- font
        'font.family': 'Times New Roman',
        'mathtext.fontset': 'cm',
        #'mathtext.fontset': 'stix',
        'font.size': 20,
        'axes.labelsize': 26,
        'legend.fontsize': 26,
        'svg.fonttype': 'path',  # Embed characters as paths
        #'svg.fonttype': 'none',  # Assume fonts are installed on the machine
        'pdf.fonttype': 42,  # embed fonts in PDF using type42 (True type)
    }
    return rcParams_dict
