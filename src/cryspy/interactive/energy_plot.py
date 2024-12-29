from ..IO import pkl_data
from ..EA.calc_hull import draw_convex_hull_binary, draw_convex_hull_ternary
from ..util.visual_util import set_params

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pymatgen.analysis.phase_diagram import PDPlotter
import numpy as np


def plot_energy(
        title=None,
        ymax=2.0,
        ymin=-0.2,
        markersize=12,
        marker_edge_width=1.0,
        marker_edge_color='black',
        alpha=1.0,
    ):
    # ---------- load data
    rslt_data = pkl_data.load_rslt()

    # ---------- EA or RS
    if 'Gen' in rslt_data.columns:
        lgen = True
    else:
        lgen = False

    # ---------- data info
    # ------ Generation
    if lgen:
        gmax = rslt_data['Gen'].max()
        print(f'Number of generation: {gmax}')
    # ------ Number of structures
    ndata = len(rslt_data)
    print(f'Number of data: {ndata}')
    # ------ check success and error
    nsuccess = rslt_data['E_eV_atom'].count()
    nerror = ndata - nsuccess
    print(f'Success: {nsuccess}')
    print(f'Error: {nerror}')
    # ------ minimum
    Emin = rslt_data['E_eV_atom'].min()
    print(f'Emin: {Emin} eV/atom')

    # ---------- setting for plot
    set_params()
    dx = 1
    xmin = -dx
    xmax = ndata + dx

    # ---------- plot
    fig, ax = plt.subplots()
    ax.set_title(title)
    ax.set_xlabel('Structure ID')
    ax.set_ylabel('Energy (eV/atom)')
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.hlines(0.0, xmin, xmax, 'k', '--')    # E = 0 line
    # ------ EA
    if lgen:
        tx = 0
        for g in range(1, gmax+1):    # generation starts from 1
            gfilter = rslt_data['Gen'] == g
            num = len(rslt_data[gfilter])
            x = np.arange(0, num) + tx
            ax.plot(
                x,
                rslt_data['E_eV_atom'][gfilter] - Emin,
                'o',
                ms=markersize,
                mew=marker_edge_width,
                markeredgecolor=marker_edge_color,
                alpha=alpha,
            )
            tx += num
    # ------ RS
    else:
        ax.plot(
            rslt_data.index,
            rslt_data['E_eV_atom'] - Emin,
            'o',
            ms=markersize,
            mew=marker_edge_width,
            markeredgecolor=marker_edge_color,
            alpha=alpha,
        )

    # ---------- return
    return fig, ax


def interact_plot_conv_hull(cgen=None, show_unstable=0.2, ternary_style='2d'):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()

    # ---------- current generation
    if cgen is None:
        cgen = max(pd_data.keys())
    pd = pd_data[cgen]

    # ---------- plot
    plotter = PDPlotter(pd, show_unstable=show_unstable, ternary_style=ternary_style)
    plotter.show()


def plot_conv_hull_binary(cgen=None, show_max=0.2, label_stable=True, vmax=0.2, bottom_margin=0.02):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()

    # ---------- current generation
    if cgen is None:
        cgen = rslt_data['Gen'].max()
    c_rslt = rslt_data[rslt_data['Gen'] == cgen]
    cgen_ids = c_rslt.index.values    # current IDs [array]
    pd = pd_data[cgen]
    hdist = hdist_data[cgen]

    # ---------- plot
    fig, ax = draw_convex_hull_binary(pd, hdist, cgen_ids, show_max, label_stable, vmax, bottom_margin)

    # ---------- return
    return fig, ax


def plot_conv_hull_ternary(cgen=None, show_max=0.2, label_stable=True, vmax=0.2):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()

    # ---------- current generation
    if cgen is None:
        cgen = rslt_data['Gen'].max()
    c_rslt = rslt_data[rslt_data['Gen'] == cgen]
    cgen_ids = c_rslt.index.values    # current IDs [array]
    pd = pd_data[cgen]
    hdist = hdist_data[cgen]

    # ---------- plot
    fig, ax = draw_convex_hull_ternary(pd, hdist, cgen_ids, show_max, label_stable, vmax)

    # ---------- return
    return fig, ax