import matplotlib.pyplot as plt

# ---------- import later
#from matplotlib.ticker import MaxNLocator
#from pymatgen.analysis.phase_diagram import PDPlotter


def set_params():
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
        'lines.markersize': 10,
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
        # ---------- save
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    }
    plt.rcParams.update(rcParams_dict)


def plot_energy_RS(rslt_data, ymax=0.2, markersize=10):
    '''
    # ---------- args
    ymax (float): max value of y-axis
    markersize (int): size of markers
    '''
    # ---------- setting
    from matplotlib.ticker import MaxNLocator
    set_params()

    # ---------- fig
    fig, ax = plt.subplots(1, 1)

    # ---------- xlim, ylim
    ndata = len(rslt_data)
    dx = max(1, ndata*0.05)    # 5% or min 1
    xmin = -dx
    xmax = ndata + dx
    ymin = -ymax*0.1
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))

    # ---------- hline at zero
    ax.axhline(y=0.0, color='k', linestyle='--')

    # ---------- plot
    ax.plot(
        rslt_data.index,
        rslt_data['E_eV_atom'] - rslt_data['E_eV_atom'].min(),
        'o',
        markersize=markersize,
        markeredgecolor='black',
    )

    # ---------- label
    ax.set_xlabel('Structure ID')
    ax.set_ylabel('Energy (eV/atom)')

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax


def plot_energy_EA(rslt_data, ymax=0.2, markersize=10):
    '''
    # ---------- args
    ymax (float): max value of y-axis
    markersize (int): size of markers
    '''
    # ---------- setting
    from matplotlib.ticker import MaxNLocator
    set_params()

    # ---------- fig
    fig, ax = plt.subplots(1, 1)

    # ---------- xlim, ylim
    ndata = len(rslt_data)
    dx = max(1, ndata*0.05)    # 5% or min 1
    xmin = -dx
    xmax = ndata + dx
    ymin = -ymax*0.1
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))

    # ---------- hline at zero
    ax.axhline(y=0.0, color='k', linestyle='--')

    # ---------- plot
    gmax = rslt_data['Gen'].max()
    for g in range(1, gmax+1):    # generation starts from 1
        gfilter = rslt_data['Gen'] == g
        ax.plot(
            rslt_data.index[gfilter],
            rslt_data['E_eV_atom'][gfilter] - rslt_data['E_eV_atom'].min(),
            'o',
            markersize=markersize,
            markeredgecolor='black',
        )

    # ---------- label
    ax.set_xlabel('Structure ID')
    ax.set_ylabel('Energy (eV/atom)')

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax


def draw_convex_hull_binary(
        phase_diagram,
        hdist,
        filtered_ids=None,
        ymax=0.2,
        label_stable=True,
        vmax=0.2,
        bottom_margin=0.04,
        markersize=10,
        axis_order='lr',
    ):
    '''
    # ---------- args
    phase_diagram (PhaseDiagram): phase diagram object
    hdist (dict): hull distance of all structures, {ID: distance, ...}
    filtered_ids (array): ID array for filtering structures to be plotted
    ymax (float): max value of y-axis
    label_stable (bool): whether to show stable compositions
    vmax (float): max value of colorbar for hull distance
    bottom_margin (float): bottom margin of y-axis
    markersize (int): size of markers
    axis_order (str): order of axis for binary phase diagram, 'lr' or 'rl'
    '''

    # ---------- setting
    set_params()
    from pymatgen.analysis.phase_diagram import PDPlotter

    # ---------- axis order
    flip_x = (axis_order == "rl")
    def tx(x):
        return 1.0 - x if flip_x else x

    # ---------- PDPlotter
    plotter_mpl = PDPlotter(
        phase_diagram,
        show_unstable=0.0,
        backend='matplotlib',
        linewidth=1.5,
        markerfacecolor='darkslateblue',
        markersize=markersize,
    )
    lines, stable_entries, unstable_entries = plotter_mpl.pd_plot_data

    # ---------- fig
    fig, ax = plt.subplots(1, 1)
    ax.set_axisbelow(True)

    # ---------- hline
    ax.axhline(y=0, xmin=0, xmax=1, color='black', linestyle='--', zorder=1)

    # ---------- draw hull lines
    for xs, ys in lines:
        xs_tx = [tx(x) for x in xs]
        ax.plot(xs_tx, ys, color="black", linewidth=1.5, zorder=1)

    # ---------- stable entries (points + optional text)
    for coord, entry in stable_entries.items():
        x = tx(coord[0])
        y = coord[1]
        ax.plot(
            x,
            y,
            'o',
            color="darkslateblue",
            markeredgecolor="black",
            markersize=markersize,
            zorder=2,
        )
        if label_stable and entry.name:
            ax.annotate(
                entry.name,
                xy=(x, y),
                xytext=(0, -10),    # 10 points vertical offset
                textcoords="offset points",
                fontsize=14,
                fontweight="normal",
                ha="center",
                va="top",
            )

    # ---------- unstable entries (colored by hull distance)
    scat_x, scat_y, scat_c = [], [], []
    s = markersize**2 / 2    # size for scatter
    for entry, coord in unstable_entries.items():
        if entry.entry_id is None:
            continue
        if filtered_ids is not None:
            if entry.entry_id not in filtered_ids:
                continue
        scat_x.append(tx(coord[0]))
        scat_y.append(coord[1])
        scat_c.append(hdist[entry.entry_id])
    mappable = ax.scatter(
        scat_x,
        scat_y,
        s=s,
        c=scat_c,
        vmin=0,
        vmax=vmax,
        cmap="Oranges_r",
        marker="D",
        edgecolors="black",
        zorder=3,
    )
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.8, pad=0.05)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Hull distance (eV/atom)', size=20, rotation=270, labelpad=30)

    # ---------- ylim
    stable_y = [coord[1] for coord in stable_entries.keys()]
    ymin = min(stable_y) - bottom_margin
    ax.set_ylim(ymin, ymax)

    # ---------- label
    # default fontweight is 'bold' in PDPlotter, so set 'normal'
    ax.set_xlabel('Composition', fontsize=20, fontweight='normal')
    ax.set_ylabel('Formation energy (eV/atom)', fontsize=20, fontweight='normal')

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax


def draw_convex_hull_ternary(
        phase_diagram,
        hdist,
        filtered_ids=None,
        show_max=0.2,
        label_stable=True,
        vmax=0.2,
        markersize=10,
        ordering=None,
    ):
    '''
    # ---------- args
    phase_diagram (PhaseDiagram): phase diagram object
    hdist (dict): hull distance of all structures, {ID: distance, ...}
    filtered_ids (array): ID array for filtering structures to be plotted
    show_max (float): max value of y-axis (binary) or hull distance (ternary)
    label_stable (bool): whether to show stable compositions
    vmax (float): max value of colorbar for hull distance
    markersize (int): size of markers
    ordering (list): order of axis for ternary phase diagram, e.g. ['A', 'B', 'C']
    '''

    # ---------- setting
    set_params()
    from pymatgen.analysis.phase_diagram import PDPlotter
    from pymatgen.analysis.phase_diagram import order_phase_diagram

    # ---------- fig
    fig, ax = plt.subplots(1, 1)
    plotter_mpl = PDPlotter(
        phase_diagram,
        show_unstable=0.0,
        backend='matplotlib',
        linewidth=1.5,
        markerfacecolor='darkslateblue',
        markersize=markersize,
    )
    plotter_mpl.get_plot(
        label_stable=label_stable,
        label_unstable=False,
        ax=ax,
        ordering=ordering,
    )

    # ---------- texts
    for text in ax.texts:
        text.set_fontsize(14)
        text.set_fontweight('normal')    # bold --> normal

    # ---------- ordering for scatter and marking
    lines, stable_entries, unstable_entries = plotter_mpl.pd_plot_data
    orderd_lines, ordered_stable_entries, ordered_unstable_entries = order_phase_diagram(
        lines, stable_entries, unstable_entries, ordering
    )

    # ---------- scatter: unstable entries
    scat_x, scat_y, scat_c = [], [], []
    s = markersize**2 / 2    # size for scatter
    for entry, coord in ordered_unstable_entries.items():
        if entry.entry_id is None:
            continue
        if filtered_ids is not None:
            if entry.entry_id not in filtered_ids:
                continue
        if hdist[entry.entry_id] <= show_max:
            scat_x.append(coord[0])
            scat_y.append(coord[1])
            scat_c.append(hdist[entry.entry_id])
    mappable = ax.scatter(
        scat_x,
        scat_y,
        s=s,
        c=scat_c,
        vmin=0,
        vmax=vmax,
        cmap='Oranges_r',
        marker='D',
        edgecolors='black',
        zorder=3,
    )
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=-0.1)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label('Hull distance (eV/atom)', size=20, rotation=270, labelpad=30)

    # ---------- return
    plt.close(fig)    # not to show the figure in Jupyter notebook when using interactive mode
    return fig, ax


def get_generation_range(plot_min_gen, plot_max_gen, hull_ref_gen, g_max_avail):
    # ------ plot_min_gen
    if plot_min_gen is not None:
        if plot_min_gen > g_max_avail:
            raise ValueError(
                f'plot_min_gen = {plot_min_gen} is larger than the maximum generation in pd_data '
                f'(latest = {g_max_avail})'
            )
        g_min = plot_min_gen
    else:
        g_min = 1

    # ------ plot_max_gen
    if plot_max_gen is not None:
        if plot_max_gen > g_max_avail:
            raise ValueError(
                f'plot_max_gen = {plot_max_gen} is larger than the maximum generation in pd_data '
                f'(latest = {g_max_avail})'
            )
        g_max = plot_max_gen
    else:
        g_max = g_max_avail

    # ------ hull_ref_gen
    if hull_ref_gen is not None:
        if hull_ref_gen > g_max_avail:
            raise ValueError(
                f'hull_ref_gen = {hull_ref_gen} is larger than the maximum generation in pd_data '
                f'Latest generation in pd_data is {g_max_avail}'
            )
        g_ref = hull_ref_gen
    elif plot_max_gen is not None:
        g_ref = plot_max_gen
    else:
        g_ref = g_max_avail

    return g_min, g_max, g_ref


def build_ordering(atype, axis_order):
    """
    Build ordering list for PDPlotter / order_phase_diagram.

    Interpretation of `axis_order`:
        Ternary (len(atype) == 3):
            axis_order is a length-3 string consisting of 't', 'l', 'r'.
            axis_order[i] indicates whether atype[i] is plotted on
            Top, Left, or Right of the ternary triangle.

    Example:
        atype = ('Li', 'Ca', 'Cl')
        axis_order = "rtl"
            Li -> Right
            Ca -> Top
            Cl -> Left
        returns ['Ca', 'Cl', 'Li']   # [Up, Left, Right]

    Returns:
        list[str] or None
            For ternary: [Up, Left, Right]
            For other dimensions: None (ordering not applied)
    """
    # ---------- ternary system
    if len(atype) == 3:
        idx_t = axis_order.index('t')
        idx_l = axis_order.index('l')
        idx_r = axis_order.index('r')
        up    = to_special_formula(atype[idx_t])
        left  = to_special_formula(atype[idx_l])
        right = to_special_formula(atype[idx_r])
        return [up, left, right]

    # ---------- higher components
    # ordering control is not supported for >3 components
    return None


def to_special_formula(sym: str) -> str:
    special_formulas = {
        "O":  "O2",
        "N":  "N2",
        "F":  "F2",
        "Cl": "Cl2",
        "H":  "H2",
    }
    return special_formulas.get(sym, sym)
