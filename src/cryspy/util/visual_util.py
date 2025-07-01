import matplotlib.pyplot as plt

# ---------- import later
#from pymatgen.analysis.phase_diagram import PDPlotter


def set_params():
    # ---------- rcParams
    rcParams_dict = {
        # ---------- figure
        'figure.figsize': (8, 6),
        'figure.dpi': 120,
        'figure.facecolor': 'white',
        # ---------- axes
        'axes.grid': True,
        'axes.linewidth': 1.5,
        'axes.labelsize': 20,
        # ---------- ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 8.0,
        'ytick.major.size': 8.0,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        # ---------- lines
        'lines.linewidth': 2.0,
        'lines.markersize': 12,
        # ---------- grid
        'grid.linestyle': ':',
        # ---------- legend
        'legend.fontsize': 20,
        # ---------- other fonts
        'font.size': 20,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Arial', 'Liberation Sans', 'DejaVu Sans', 'sans'],
        'mathtext.fontset': 'cm',
        #'mathtext.fontset': 'stix',
        'svg.fonttype': 'path',  # Embed characters as paths
        #'svg.fonttype': 'none',  # Assume fonts are installed on the machine
        'pdf.fonttype': 42,  # embed fonts in PDF using type42 (True type)
    }
    plt.rcParams.update(rcParams_dict)


def draw_convex_hull_binary(
        pd,
        hdist,
        cgen_ids,
        show_max=0.2,
        label_stable=True,
        vmax=0.2,
        bottom_margin=0.02,
    ):
    '''
    # ---------- args
    pd (PhaseDiagram): phase diagram object
    hdist (dict): hull distance of all structures, {ID: distance, ...}
    cgen_ids (array): ID array of current generation structures
    show_max (float): max value of y-axis (binary) or hull distance (ternary)
    label_stable (bool): whether to show stable compositions
    vmax (float): max value of colorbar for hull distance
    bottom_margin (float): bottom margin of y-axis
    '''

    # ---------- setting
    set_params()
    from pymatgen.analysis.phase_diagram import PDPlotter

    # ---------- fig
    fig, ax = plt.subplots(1, 1)
    plotter_mpl = PDPlotter(pd, show_unstable=0.0, backend='matplotlib', linewidth=1.5, markerfacecolor='darkslateblue', markersize=10)
    plotter_mpl.get_plot(label_stable=label_stable, label_unstable=False, ax=ax)
    ax.set_axisbelow(True)

    # ---------- hline
    ax.axhline(y=0, xmin=0, xmax=1, color='black', linestyle='--', zorder=1)

    # ---------- label for only binary system
    # default fontweight is 'bold' in PDPlotter, so set 'normal'
    ax.set_xlabel('Composition', fontsize=20, fontweight='normal')
    ax.set_ylabel('Formation energy (eV/atom)', fontsize=20, fontweight='normal')

    # ---------- texts
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_fontweight('normal')    # bold --> normal

    # ---------- scatter: unstable entries
    scat_x = []
    scat_y= []
    scat_c = []
    lines, stable_entries, unstable_entries = plotter_mpl.pd_plot_data
    for entry, coord in unstable_entries.items():
        if entry.entry_id is not None:
            scat_x.append(coord[0])
            scat_y.append(coord[1])
            scat_c.append(hdist[entry.entry_id])
    mappable = ax.scatter(scat_x, scat_y, s=50, c=scat_c, vmin=0, vmax=vmax, cmap='Oranges_r', marker='D', edgecolors='black', zorder=2)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, pad=0.05)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Hull distance (eV/atom)', size=20, rotation=270, labelpad=30)

    # ---------- mark the current generation
    stable_compos = {entry.entry_id: compos for compos, entry in stable_entries.items()}
    unstable_compos = {entry.entry_id: compos for entry, compos in unstable_entries.items()}
    for cid in cgen_ids:
        if cid in stable_compos:
            mx, my = stable_compos[cid][0], stable_compos[cid][1]
            ax.plot(mx, my, '+', markeredgecolor='white')
        elif cid in unstable_compos:
            mx, my = unstable_compos[cid][0], unstable_compos[cid][1]
            ax.plot(mx, my, '+', markersize=10, markeredgewidth=0.5,  markeredgecolor='navy')

    # ---------- ylim
    stable_y = list(stable_entries.keys())
    ymin = min(stable_y, key=lambda x: x[1])[1] - bottom_margin
    ax.set_ylim(ymin, show_max)

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax


def draw_convex_hull_ternary(
        pd,
        hdist,
        cgen_ids,
        show_max=0.2,
        label_stable=True,
        vmax=0.2,
    ):
    '''
    # ---------- args
    pd (PhaseDiagram): phase diagram object
    hdist (dict): hull distance of all structures, {ID: distance, ...}
    cgen_ids (array): ID array of current generation structures
    show_max (float): max value of y-axis (binary) or hull distance (ternary)
    label_stable (bool): whether to show stable compositions
    vmax (float): max value of colorbar for hull distance
    '''

    # ---------- setting
    set_params()
    from pymatgen.analysis.phase_diagram import PDPlotter

    # ---------- fig
    fig, ax = plt.subplots(1, 1)
    plotter_mpl = PDPlotter(pd, show_unstable=0.0, backend='matplotlib', linewidth=1.5, markerfacecolor='darkslateblue', markersize=10)
    plotter_mpl.get_plot(label_stable=label_stable, label_unstable=False, ax=ax)

    # ---------- texts
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_fontweight('normal')    # bold --> normal

    # ---------- scatter: unstable entries
    scat_x = []
    scat_y= []
    scat_c = []
    lines, stable_entries, unstable_entries = plotter_mpl.pd_plot_data
    for entry, coord in unstable_entries.items():
        if entry.entry_id is not None:
            if hdist[entry.entry_id] <= show_max:
                scat_x.append(coord[0])
                scat_y.append(coord[1])
                scat_c.append(hdist[entry.entry_id])
    mappable = ax.scatter(scat_x, scat_y, s=30, c=scat_c, vmin=0, vmax=vmax, cmap='Oranges_r', marker='D', edgecolors='black', zorder=3)
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=-0.1)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Hull distance (eV/atom)', size=20, rotation=270, labelpad=30)

    # ---------- mark the current generation
    stable_compos = {entry.entry_id: compos for compos, entry in stable_entries.items()}
    unstable_compos = {entry.entry_id: compos for entry, compos in unstable_entries.items()}
    for cid in cgen_ids:
        if cid in hdist and hdist[cid] <= show_max:
            if cid in stable_compos:
                mx, my = stable_compos[cid][0], stable_compos[cid][1]
                ax.plot(mx, my, '+', markeredgecolor='white', zorder=2)
            elif cid in unstable_compos:
                mx, my = unstable_compos[cid][0], unstable_compos[cid][1]
                ax.plot(mx, my, '+', markersize=6, markeredgewidth=0.5,  markeredgecolor='navy', zorder=4)

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax
