from logging import getLogger

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import ConvexHull


logger = getLogger('cryspy')


def calc_convex_hull(atype, ratio_data, ef_all, c_ids, gen, vmax=None):
    '''
    Input:
        atype (tuple): atom type, e.g. ('Na', 'Cl')
        ratio_data [dict]: ratio of all structures, {ID: (ratio tuple), ...}
        ef_all [dict]: formation energy of all structures, {ID: Ef, ...}
        c_ids [array]: ID array of current generation structures, just for plot
        gen [int]: current generation, just for plot
        vmax [float]: max value of colorbar for hull distance, just for plot

    Return:
        hdist [dict]: hull distance of all structures, {ID: distance, ...}
    '''

    # ---------- remove Ef >= 0.0
    ratio_chull = []
    ef_chull = []
    for cid in ef_all:
        if ef_all[cid] < 0.0:
            ratio_chull.append(ratio_data[cid])
            ef_chull.append(ef_all[cid])

    # ---------- no negative Ef
    #            in this case, hull distance is equivalent to Ef
    if not ef_chull:
        # np.nan --> np.inf in hdist
        hdist = {cid: np.inf if np.isnan(ef_all[cid]) else ef_all[cid] for cid in ef_all}
        if len(atype) == 2:
            draw_convex_hull_2d(atype, ratio_data, ef_all, c_ids, gen, None)
        elif len(atype) == 3:
            draw_convex_hull_3d(atype, ratio_data, hdist, c_ids, gen, None, vmax)
        return hdist

    # ---------- add end points
    n = len(atype)
    end_ratio = [[1 if i == j else 0 for j in range(n)] for i in range(n)]    # identity_matrix
    end_ef = [0.0 for i in range(n)]
    # ------ extend
    ratio_chull.extend(end_ratio)    # last n data are end points
    ef_chull.extend(end_ef)

    # ---------- to array
    ratio_chull = np.array(ratio_chull)
    ratio_chull = ratio_chull[:, :-1]    # last composition is not used
    ef_chull = np.array(ef_chull)
    points = np.hstack([ratio_chull, ef_chull[:, np.newaxis]])

    # ---------- calc convex hull
    hull = ConvexHull(points)
    vpoints = hull.points[hull.vertices]    # hull.vertices: index of vertices in hull.points

    # ---------- hull distance
    # ------ binary
    if len(atype) == 2:
        sorted_indices = np.argsort(vpoints[:, 0])
        sorted_vx = vpoints[:, 0][sorted_indices]
        sorted_vef = vpoints[:, 1][sorted_indices]
        x = [ratio_data[cid][0] for cid in ef_all]    # first composition ratio, order: ID in ef_all
        ef_on_hull = np.interp(x, sorted_vx, sorted_vef)    # order: ID in ef_all
    # ------ ternary or more
    else:
        interp = LinearNDInterpolator(vpoints[:, :-1], vpoints[:, -1])    # no need to sort
        x = [(ratio_data[cid][:-1]) for cid in ef_all]    # (n-1) dim. composition ratio, order: ID in ef_all
        ef_on_hull = interp(x)    # order: ID in ef_all
    # ------ calc hull distance
    hdist = {
        cid: np.inf if np.isnan(ef_all[cid])
        else ef_all[cid] - ef_on_hull[i]
        for i, cid in enumerate(ef_all)
    }

    # ---------- draw convex hull
    if len(atype) == 2:
        draw_convex_hull_2d(atype, ratio_data, ef_all, c_ids, gen, hull)
    elif len(atype) == 3:
        draw_convex_hull_3d(atype, ratio_data, hdist, c_ids, gen, hull, vmax)

    # ---------- return
    return hdist


def draw_convex_hull_2d(atype, ratio_data, ef_all, c_ids, gen, hull):
    '''
    # ---------- args
    atype (tuple): atom type, e.g. ('Na', 'Cl')
    ratio_data [dict]: ratio of all structures, {ID: (ratio tuple), ...}
    ef_all [dict]: formation energy of all structures, {ID: Ef, ...}
    c_ids [array]: ID array of current generation structures
    gen [int]: current generation
    hull [ConvexHull]: convex hull object
    '''

    # ---------- setting
    plt.rcParams.update(_set_params())

    # ---------- fig
    fig, ax = plt.subplots(1, 1)

    # ---------- hline
    ax.axhline(y=0, xmin=0, xmax=1, color='black', linestyle='--')

    # ---------- label
    ax.set_xlabel('$x$ in '+f'{atype[0]}'+'$_{x}$'+f'{atype[1]}'+'$_{1-x}$')
    ax.set_ylabel('Formation energy (eV/atom)')

    # ---------- lim
    min_ef = min(ef_all.values())
    ymin = min_ef - 0.01 if min_ef < 0 else -0.01
    ax.set_ylim(ymin, 0.05)

    # ---------- plot ef_all
    for cid in ef_all:
        ax.plot(ratio_data[cid][0], ef_all[cid], 'o', color='thistle', markeredgecolor='navy', alpha=0.75)

    # ---------- plot convex hull
    if hull is None:    # no negative Ef, plot only end points
        ax.plot([0, 1], [0, 0], 'o', color='darkslateblue', markeredgecolor='navy')
    else:
        for i, simplex in enumerate(hull.simplices):
            if i != 0:    # i == 0 --> remove line from (1, 0) to (0, 0)
                ax.plot(hull.points[simplex, 0], hull.points[simplex, 1], '-o', color='darkslateblue', markeredgecolor='navy')

    # ---------- plot current generation
    for cid in c_ids:
        ax.plot(ratio_data[cid][0], ef_all[cid], '+', markeredgecolor='hotpink')

    # ---------- save figure
    fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.png', bbox_inches='tight')


def draw_convex_hull_3d(atype, ratio_data, hdist, c_ids, gen, hull, vmax=None):
    '''
    # ---------- args
    atype (tuple): atom type, e.g. ('Na', 'Cl')
    ratio_data [dict]: ratio of all structures, {ID: (ratio tuple), ...}
    hdist [dict]: hull distance of all structures, {ID: distance, ...}
    c_ids [array]: ID array of current generation structures
    gen [int]: current generation
    hull [ConvexHull]: convex hull object
    vmax [float]: max value of colorbar for hull distance
    '''

    # ---------- setting
    plt.rcParams.update(_set_params())

    # ---------- fig
    fig, ax = plt.subplots(1, 1)
    ax.axis('equal')
    ticks = np.arange(0, 1.1, 0.2)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    # ---------- label
    ax.set_xlabel('$x$ in ' + f'{atype[0]}'+'$_{x}$' + f'{atype[1]}'+'$_{y}$' + f'{atype[2]}'+'$_{z}$')
    ax.set_ylabel('$y$ in ' + f'{atype[0]}'+'$_{x}$' + f'{atype[1]}'+'$_{y}$' + f'{atype[2]}'+'$_{z}$')

    # ---------- plot convex hull
    if hull is None:    # no negative Ef, plot only end points
        ax.plot([0, 0, 1, 0], [0, 1, 0, 0], '-o', color='darkslateblue', markeredgecolor='navy')
    else:
        for i, simplex in enumerate(hull.simplices):
            spoints_closed = np.vstack([hull.points[simplex], hull.points[simplex[0]]])    # append the first point to close the line
            ax.plot(spoints_closed[:, 0], spoints_closed[:, 1], '-o', color='darkslateblue', markeredgecolor='navy')

    # ---------- hdist
    #            scatter is drawn on the back of the normal plot
    vmin = 0
    cmap = 'Purples_r'
    x = [ratio_data[cid][0] for cid in hdist]
    y = [ratio_data[cid][1] for cid in hdist]
    c = list(hdist.values())
    mappable = ax.scatter(x, y, c=c, vmin=vmin, vmax=vmax, cmap=cmap, edgecolors='navy', alpha=0.75)
    fig.colorbar(mappable, ax=ax)

    # ---------- current generation
    for cid in c_ids:
        ax.plot(ratio_data[cid][0], ratio_data[cid][1], '+', markeredgecolor='hotpink')

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
        'lines.linewidth': 1.5,
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










# def calc_convex_hull_2d(rin, ratio_data, ef_all, c_ids, gen):
#     '''
#     Input:
#         ratio_data [dict]: ratio of all structures, {ID: (ratio tuple), ...}
#         ef_all [dict]: formation energy of all structures, {ID: Ef, ...}
#         c_ids [array]: ID array of current generation structures
#         gen [int]: current generation

#     Return:
#         hdist [dict]: hull distance of all structures, {ID: distance, ...}
#     '''
#     # ---------- initialize
#     ratio_chull =  [0.0, 1.0]    # only ef < 0 for calculation of convex hull
#     ef_chull =  [0.0, 0.0]       # only ef < 0 for calculation of convex hull

#     # ---------- dict to list
#     for cid in ratio_data:
#         if ef_all[cid] < 0.0:
#             ratio_chull.append(ratio_data[cid][0])
#             ef_chull.append(ef_all[cid])

#     # ---------- no negative Ef
#     #            in this case, hull distance is equivalent to Ef
#     if len(ef_chull) == 2:    # only end points
#         # np.nan --> np.inf in hdist
#         hdist = {cid: np.inf if np.isnan(ef_all[cid]) else ef_all[cid] for cid in ef_all}
#         draw_convex_hull_2d(rin, None, ratio_data, ef_all, c_ids, gen)
#         return hdist

#     # ---------- calc convex hull
#     points = list(zip(ratio_chull, ef_chull))
#     hull = ConvexHull(points)
#     vpoints = hull.points[hull.vertices]    # hull.vertices: index of vertices in hull.points
#     vpoints = np.vstack((vpoints, vpoints[0]))    # just for plot. vpoints[0] should be [1, 0] by ConvexHull()

#     # ---------- hull distance
#     hdist = {
#         cid: np.inf if np.isnan(ef_all[cid])
#         else hull_distance_2d(ratio_data[cid][0], ef_all[cid], hull.equations)
#         for cid in ef_all
#     }

#     # ---------- draw convex hull
#     draw_convex_hull_2d(rin, vpoints, ratio_data, ef_all, c_ids, gen)

#     # ---------- return
#     return hdist



# def hull_distance_2d(x, y, equations):
#     '''
#     equations: [eq0, eq1, eq2, ...], eq0: [a, b, c] for a*x + b*y + c = 0

#     Find the distance from all equations and adopt the minimum distance.
#     '''
#     hdists = []
#     for eq in equations:
#         if eq[0] != 0.0:
#             dist = y - (-eq[0]*x - eq[2])/eq[1]
#             if dist < 0.0:
#                 logger.warning('hdist <= 0.0, check hull distance.')
#             hdists.append(dist)

#     return min(hdists)