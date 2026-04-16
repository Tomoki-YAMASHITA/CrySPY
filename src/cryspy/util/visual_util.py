import matplotlib.pyplot as plt
import numpy as np
import os

# ---------- import later
#from matplotlib.ticker import MaxNLocator
#import mpltern
#from pymatgen.analysis.phase_diagram import PDPlotter

#from .struc_util import get_feasible_composition


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
    max_indx = rslt_data.index.max()
    dx = max(1, max_indx*0.05)    # 5% or min 1
    xmin = -dx
    xmax = max_indx + dx
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


def plot_energy_EA(rslt_data, ymax=0.2, markersize=10, max_indx=None):
    '''
    # ---------- args
    ymax (float): max value of y-axis
    markersize (int): size of markers
    max_indx (int): max value of index for x-axis
    '''
    # ---------- setting
    from matplotlib.ticker import MaxNLocator
    set_params()

    # ---------- fig
    fig, ax = plt.subplots(1, 1)

    # ---------- xlim, ylim
    if max_indx is None:
        max_indx = rslt_data.index.max()
    dx = max(1, max_indx*0.05)    # 5% or min 1
    xmin = -dx
    xmax = max_indx + dx
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
        show_max=None,
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
    show_max (float): Maximum hull distance for entries to be shown in the plot.
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
        if show_max is None or hdist[entry.entry_id] <= show_max:
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


def get_generation_range(plot_min_gen, plot_max_gen, ref_gen, g_max_avail):
    # ------ plot_min_gen
    if plot_min_gen is not None:
        if plot_min_gen > g_max_avail:
            raise ValueError(
                f'plot_min_gen = {plot_min_gen} is larger than the latest generation '
                f'(latest = {g_max_avail})'
            )
        g_min = plot_min_gen
    else:
        g_min = 1

    # ------ plot_max_gen
    if plot_max_gen is not None:
        if plot_max_gen > g_max_avail:
            raise ValueError(
                f'plot_max_gen = {plot_max_gen} is larger than the latest generation '
                f'(latest = {g_max_avail})'
            )
        g_max = plot_max_gen
    else:
        g_max = g_max_avail

    # ------ ref_gen
    if ref_gen is not None:
        if ref_gen > g_max_avail:
            raise ValueError(
                f'ref_gen = {ref_gen} is larger than the latest generation '
                f'(latest = {g_max_avail})'
            )
        g_ref = ref_gen
    else:
        g_ref = g_max_avail

    return g_min, g_max, g_ref


def build_ordering(atype, axis_order, use_special_formula=True):
    """
    Build ordering list for PDPlotter / order_phase_diagram.

    # ---------- Interpretation of `axis_order`:
        Ternary (len(atype) == 3):
            axis_order is a length-3 string consisting of 't', 'l', 'r'.
            axis_order[i] indicates whether atype[i] is plotted on
            Top, Left, or Right of the ternary triangle.

    # ---------- Example:
        atype = ('Li', 'Ca', 'Cl')
        axis_order = "rtl"
            Li -> Right
            Ca -> Top
            Cl -> Left
        returns ['Ca', 'Cl', 'Li']   # [Up, Left, Right]

    # ---------- Parameters
        atype : sequence[str]
            Component symbols, e.g. ('Li', 'Ca', 'Cl').
        axis_order : str
            Permutation of 't', 'l', 'r' (top / left / right).
        use_special_formula : bool, optional
            If True, apply to_special_formula (e.g. O -> O2, Cl -> Cl2).
            If False, use raw symbols in atype without conversion.

    # ---------- Returns:
        list[str] or None
            For ternary: [Up, Left, Right]
            For other dimensions: None (ordering not applied)
    """
    # ---------- ternary system
    if len(atype) == 3:
        idx_t = axis_order.index('t')
        idx_l = axis_order.index('l')
        idx_r = axis_order.index('r')
        if use_special_formula:
            up    = to_special_formula(atype[idx_t])
            left  = to_special_formula(atype[idx_l])
            right = to_special_formula(atype[idx_r])
        else:
            up    = atype[idx_t]
            left  = atype[idx_l]
            right = atype[idx_r]
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


def plot_composition_window_2d(
        atype,
        min_comp,
        max_comp,
        tol=1e-12,
        ndigits=2,
        axis_order='lr',
    ):
    """
    Visualize composition constraints for a binary system.

    Always plots the two input composition bands.
    If an overlap (feasible composition range) exists, its bounds are shown
    as vertical dashed lines with numeric labels.

    Parameters
    ----------
    atype : tuple[str, str]
        Element names, e.g. ('Cu', 'Au').
    min_comp : tuple[float, float]
        Minimum composition fractions for each component.
    max_comp : tuple[float, float]
        Maximum composition fractions for each component.
    tol : float, optional
        Numerical tolerance used when judging feasibility and overlap.
    ndigits : int, optional
        Number of decimal digits used when displaying overlap bounds
        in the figure (display only; does not affect calculations).
    axis_order : str, optional
        Order of axis for binary system, 'lr' or 'rl'.
        'lr': atype[0] is on the left and atype[1] is on the right.
        'rl': atype[1] is on the left and atype[0] is on the right.
    """

    # ---------- setting
    set_params()
    from .struc_util import get_feasible_composition

    # ---------- axis order
    flip_x = (axis_order == 'rl')

    def tx(x):
        return 1.0 - x if flip_x else x

    def tx_interval(xmin, xmax):
        return sorted((tx(xmin), tx(xmax)))

    # ---------- visual parameters (internal constants)
    alpha_base = 0.6
    alpha_overlap = 1.0

    # ---------- extract component bounds
    x1_min = float(min_comp[0])
    x1_max = float(max_comp[0])
    x1_from_c2_min = 1.0 - float(max_comp[1])
    x1_from_c2_max = 1.0 - float(min_comp[1])

    # ---------- transform intervals for display
    x1_plot_min, x1_plot_max = tx_interval(x1_min, x1_max)
    x2_plot_min, x2_plot_max = tx_interval(x1_from_c2_min, x1_from_c2_max)

    # ---------- figure setup
    fig, ax = plt.subplots(figsize=(8, 1.5))

    # ---------- base bands (always)
    ax.fill_between([x1_plot_min, x1_plot_max], -0.5, 0.0,
                    color="C0", alpha=alpha_base)
    ax.fill_between([x2_plot_min, x2_plot_max], 0.0, 0.5,
                    color="C1", alpha=alpha_base)

    # ---------- element labels
    ax.text(tx(x1_min + 0.05), -0.25, atype[0], ha="center", va="center")
    ax.text(tx(x1_from_c2_max - 0.05), 0.25, atype[1], ha="center", va="center")

    # ---------- feasible composition range
    feasible_comp = get_feasible_composition(min_comp, max_comp, tol=tol)
    if feasible_comp is not None:
        overlap_min, overlap_max = feasible_comp[0]
        overlap_plot_min, overlap_plot_max = tx_interval(overlap_min, overlap_max)

        ax.fill_between([overlap_plot_min, overlap_plot_max], -0.5, 0.0,
                        color="C0", alpha=alpha_overlap)
        ax.fill_between([overlap_plot_min, overlap_plot_max], 0.0, 0.5,
                        color="C1", alpha=alpha_overlap)
        ax.axvline(overlap_plot_min, color="k", linestyle="--", linewidth=1)
        ax.axvline(overlap_plot_max, color="k", linestyle="--", linewidth=1)
        ax.text(overlap_plot_min, 0.55, f"{overlap_plot_min:.{ndigits}f}",
                ha="center", va="bottom")
        ax.text(overlap_plot_max, 0.55, f"{overlap_plot_max:.{ndigits}f}",
                ha="center", va="bottom")
    else:
        ax.text(0.5, 0.55, "Infeasible", ha="center", va="bottom")

    # ---------- axis settings
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Composition")

    # ---------- return
    return fig, ax


def plot_composition_window_3d(
        atype,
        min_comp,
        max_comp,
        axis_order='tlr',
        tol=1e-12,
    ):
    """
    Visualize composition constraints for a ternary system.

    Always plots the feasible region on the ternary diagram.
    If an overlap (feasible composition region) exists,
    it is highlighted with color.

    Parameters
    ----------
    atype : tuple[str, str, str]
        Element names, e.g. ('Cu', 'Au', 'Ag').
    min_comp : tuple[float, float, float]
        Minimum composition fractions for each component.
    max_comp : tuple[float, float, float]
        Maximum composition fractions for each component.
    axis_order : str, optional
        Permutation of 't', 'l', 'r' specifying which component
        goes to top / left / right. Default is 'tlr'.
        For example, 'rtl' means:
            atype[0] -> right
            atype[1] -> top
            atype[2] -> left
    tol : float, optional
        Numerical tolerance used when judging feasibility and overlap.
    """

    # ---------- setting
    set_params()
    import mpltern

    # ---------- axis order
    idx_t = axis_order.index('t')
    idx_l = axis_order.index('l')
    idx_r = axis_order.index('r')

    def to_tlr(x0, x1, x2):
        comps = (x0, x1, x2)
        return comps[idx_t], comps[idx_l], comps[idx_r]

    # ---------- figure setup
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='ternary')

    # ---------- draw bands for each component
    for i in range(3):
        for c in (min_comp[i], max_comp[i]):
            line = _get_feasible_segment(i, c, tol=tol)
            if line is None:
                continue
            t, l, r = to_tlr(*line)
            ax.plot(t, l, r, color='k', linestyle="--", linewidth=1.0)

    # ---------- feasible polygon
    poly = _get_feasible_polygon_ternary(min_comp, max_comp, tol=tol)
    if poly is None:
        ax.text(1/3, 1/3, 1/3, "Infeasible", ha="center", va="center")
    else:
        t, l, r = to_tlr(*poly)
        n = len(t)
        if n >= 3:
            # ------ 2D feasible polygon
            ax.fill(
                t, l, r,
                color="C0",
                alpha=0.6,
                edgecolor="k",
            )
        elif n == 2:
            # ------ 1D feasible line segment
            ax.plot(
                t, l, r,
                color="C0",
                linewidth=3.0,
            )
        elif n == 1:
            # ------ 0D feasible single composition
            ax.plot(
                t, l, r,
                marker="o",
                markersize=8,
                color="C0",
                markeredgecolor="k",
            )

    # ---------- axis labels (mpltern)
    ax.set_tlabel(atype[idx_t])
    ax.set_llabel(atype[idx_l])
    ax.set_rlabel(atype[idx_r])

    # ---------- return
    return fig, ax


def _get_feasible_segment(i, c, tol=1e-12):
    """
    Helper function to compute the line segment corresponding to
    a constant composition constraint x_i = c on a ternary simplex.

    We consider a ternary composition expressed as:
        x1 + x2 + x3 = 1
        x1 >= 0, x2 >= 0, x3 >= 0

    This function returns the two end points of the line segment
    where one component (xi) is fixed to a constant value.

    Parameters
    ----------
    i : int
        Component index (0, 1, or 2) corresponding to x1, x2, x3.
    c : float
        Fixed composition value (assumed 0 <= c <= 1).
    tol : float
        Numerical tolerance.

    Returns
    -------
    x1, x2, x3 : np.ndarray or None
        Each is a 1D array of length 2 representing the two end points
        of the line segment on the simplex.
        If the line is effectively outside the simplex
        (e.g. c <= 0 or c >= 1), returns None.
    """
    c = float(c)
    if c <= 0.0 + tol or c >= 1.0 - tol:
        return None

    if i == 0:
        # x1 = c, x2 + x3 = 1 - c
        x1 = np.array([c, c])
        x2 = np.array([0.0, 1.0 - c])
        x3 = 1.0 - x1 - x2

    elif i == 1:
        # x2 = c, x1 + x3 = 1 - c
        x2 = np.array([c, c])
        x1 = np.array([0.0, 1.0 - c])
        x3 = 1.0 - x1 - x2

    elif i == 2:
        # x3 = c, x1 + x2 = 1 - c
        x3 = np.array([c, c])
        x1 = np.array([0.0, 1.0 - c])
        x2 = 1.0 - x1 - x3

    else:
        raise ValueError("i must be 0, 1, or 2")

    # Safety: clip within [0, 1] to avoid tiny numerical overflow
    x1 = np.clip(x1, 0.0, 1.0)
    x2 = np.clip(x2, 0.0, 1.0)
    x3 = np.clip(x3, 0.0, 1.0)

    return x1, x2, x3


def _get_feasible_polygon_ternary(min_comp, max_comp, tol=1e-12):
    """
    Compute the feasible composition region for a ternary system under
    lower/upper composition bounds, and return its polygon vertices.

    This function assumes a ternary composition expressed as:
        x1 + x2 + x3 = 1
    where
        x1 : composition of component 1
        x2 : composition of component 2
        x3 : composition of component 3 (= 1 - x1 - x2)

    Given
        min_comp[i] <= xi <= max_comp[i]    (i = 1, 2, 3)
    the intersection of these constraints with the ternary simplex
    (x1 >= 0, x2 >= 0, x3 >= 0) is computed as a polygon.

    Parameters
    ----------
    min_comp : array-like of length 3
        Lower bounds of compositions (x1_min, x2_min, x3_min)
    max_comp : array-like of length 3
        Upper bounds of compositions (x1_max, x2_max, x3_max)
    tol : float, optional
        Numerical tolerance used when judging feasibility and intersection.

    Returns
    -------
    None
        If no feasible composition region exists.
    (x1, x2, x3) : tuple of np.ndarray
        Vertex coordinates of the feasible polygon.
        Each of x1, x2, x3 is a 1D array of shape (n_vertices,).
        Vertices are sorted in counter-clockwise order.
    """
    # ---------- import
    from .struc_util import get_feasible_composition

    # ---------- convert to numpy arrays for easier calculations
    min_comp = np.array(min_comp, float)
    max_comp = np.array(max_comp, float)

    # ---------- check feasibility
    feasible = get_feasible_composition(min_comp, max_comp, tol=tol)
    if feasible is None:
        return None

    # ---------- Outer boundary of the ternary composition triangle
    lines = []
    lines.append(("x1", 0.0))       # x1 = 0
    lines.append(("x2", 0.0))       # x2 = 0
    lines.append(("x1+x2", 1.0))    # x3 = 0 --> x1 + x2 = 1

    # ---------- Feasible band (lower/upper bounds inside the triangle)
    # ------ x1
    if min_comp[0] > 0.0 + tol:
        lines.append(("x1", float(min_comp[0])))
    if max_comp[0] < 1.0 - tol:
        lines.append(("x1", float(max_comp[0])))
    # ------ x2
    if min_comp[1] > 0.0 + tol:
        lines.append(("x2", float(min_comp[1])))
    if max_comp[1] < 1.0 - tol:
        lines.append(("x2", float(max_comp[1])))
    # ----- x3 --> 1.0 - (x1 + x2)
    if min_comp[2] > 0.0 + tol:
        c = 1.0 - float(min_comp[2])   # x1 + x2 <= c
        lines.append(("x1+x2", c))
    if max_comp[2] < 1.0 - tol:
        c = 1.0 - float(max_comp[2])   # x1 + x2 >= c
        lines.append(("x1+x2", c))

    # ---------- Compute all candidate intersection points of boundary lines
    pts = []
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            # ------ intersection of two constraint lines
            p = _intersect_two_lines(lines[i], lines[j], tol=tol)
            if p is None:
                continue
            x1, x2 = p
            x3 = 1.0 - x1 - x2
            # ------ Basic simplex validity check (inside triangle)
            #        must satisfy x_i >= 0 and x_i <= 1
            if x1 < -tol or x2 < -tol or x3 < -tol:
                continue
            if x1 > 1.0 + tol or x2 > 1.0 + tol or x3 > 1.0 + tol:
                continue
            # ------ Check min/max composition constraints
            if not (min_comp[0] - tol <= x1 <= max_comp[0] + tol):
                continue
            if not (min_comp[1] - tol <= x2 <= max_comp[1] + tol):
                continue
            if not (min_comp[2] - tol <= x3 <= max_comp[2] + tol):
                continue
            # ------ this intersection point is feasible --> keep it
            pts.append((x1, x2))
    # ------ if no feasible vertices exist, feasible region does not exist
    if not pts:
        return None

    # ---------- Remove duplicate points (caused by numerical tolerance)
    uniq = []
    for x1, x2 in pts:
        if all((abs(x1 - u1) > tol or abs(x2 - u2) > tol) for (u1, u2) in uniq):
            uniq.append((x1, x2))

    # ---------- Sort points in counter-clockwise order
    pts_arr = np.array(uniq)
    cx, cy = pts_arr.mean(axis=0)    # centroid
    angles = np.arctan2(pts_arr[:, 1] - cy, pts_arr[:, 0] - cx)
    order = np.argsort(angles)
    pts_sorted = pts_arr[order]

    # ---------- sorted points
    x1 = pts_sorted[:, 0]
    x2 = pts_sorted[:, 1]
    x3 = 1.0 - x1 - x2

    # ---------- return
    return x1, x2, x3



def _intersect_two_lines(l1, l2, tol=1e-12):
    """
    Compute the intersection of two constraint lines in 2D composition space.

    Each line is represented in one of the following forms:
        ("x1", c)       →  x1 = c
        ("x2", c)       →  x2 = c
        ("x1+x2", c)    →  x1 + x2 = c   (equivalently x3 = 1 - c in ternary)

    Parameters
    ----------
    l1, l2 : tuple
        Line specification in the above format.
    tol : float
        Numerical tolerance (not heavily used here but kept for consistency).

    Returns
    -------
    (x1, x2) : tuple[float, float]
        Intersection point in (x1, x2) coordinates.
    None :
        If the two lines are parallel or otherwise do not intersect
        in a meaningful way.
    """
    t1, c1 = l1
    t2, c2 = l2

    # ---------- Case 1: Same type of constraint --> parallel, no intersection
    if t1 == t2:
        return None

    # ---------- Case 2: x1 = const  AND  x2 = const
    #            Directly determined without solving equations
    if {"x1", "x2"} == {t1, t2}:
        x1 = c1 if t1 == "x1" else c2
        x2 = c2 if t2 == "x2" else c1
        return (x1, x2)

    # ---------- Case 3: x1 = const  AND  (x1 + x2 = const)
    #            substitute x1 into the sum constraint
    if (t1 == "x1" and t2 == "x1+x2") or (t2 == "x1" and t1 == "x1+x2"):
        x1 = c1 if t1 == "x1" else c2
        s  = c2 if t2 == "x1+x2" else c1
        x2 = s - x1
        return (x1, x2)

    # ---------- Case 4: x2 = const  AND  (x1 + x2 = const)
    #            substitute x2 into the sum constraint
    if (t1 == "x2" and t2 == "x1+x2") or (t2 == "x2" and t1 == "x1+x2"):
        x2 = c1 if t1 == "x2" else c2
        s  = c2 if t2 == "x1+x2" else c1
        x1 = s - x2
        return (x1, x2)

    # ---------- This should not happen under current design
    return None


def save_composition_window(
        atype,
        gen,
        min_comp,
        max_comp,
        fig_format,
        axis_order=None,
    ):
    if min_comp is None and max_comp is None:
        return

    if len(atype) == 2:
        if axis_order is None:
            axis_order = 'lr'
        fig, _ = plot_composition_window_2d(
            atype=atype,
            min_comp=min_comp,
            max_comp=max_comp,
            axis_order=axis_order,
        )
    elif len(atype) == 3:
        if axis_order is None:
            axis_order = 'tlr'
        fig, _ = plot_composition_window_3d(
            atype=atype,
            min_comp=min_comp,
            max_comp=max_comp,
            axis_order=axis_order,
        )
    else:
        return

    os.makedirs('./data/convex_hull', exist_ok=True)
    fig.savefig(f'./data/convex_hull/composition_window_{gen}.{fig_format}')
