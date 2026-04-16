from ..IO import pkl_data
from ..util.visual_util  import (
    plot_energy_RS,
    plot_energy_EA,
    draw_convex_hull_binary,
    draw_convex_hull_ternary,
    get_generation_range,
    build_ordering
)
from pymatgen.analysis.phase_diagram import PDPlotter


def plot_energy(ymax=0.2, markersize=10):
    # ---------- load data
    rslt_data = pkl_data.load_rslt()

    # ---------- EA or RS
    if 'Gen' in rslt_data.columns:
        lgen = True
    else:
        lgen = False

    # ---------- plot
    if not lgen:    # RS
        fig, ax = plot_energy_RS(rslt_data, ymax, markersize)
    else:    # EA
        fig, ax = plot_energy_EA(rslt_data, ymax, markersize)

    # ---------- return
    return fig, ax


def interactive_plot_convex_hull(cgen=None, show_unstable=0.2, ternary_style='2d'):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()

    # ---------- current generation
    if cgen is None:
        cgen = max(pd_data.keys())
    phase_diagram = pd_data[cgen]

    # ---------- plot
    plotter = PDPlotter(phase_diagram, show_unstable=show_unstable, ternary_style=ternary_style)
    plotter.show()


def plot_convex_hull_binary(
    plot_min_gen=None,
    plot_max_gen=None,
    ref_gen=None,
    ymax=0.2,
    label_stable=True,
    vmax=0.2,
    bottom_margin=0.04,
    markersize=10,
    axis_order='lr',
):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()

    # ---------- generation range
    g_max_avail = max(pd_data.keys())
    g_min, g_max, g_ref = get_generation_range(
        plot_min_gen=plot_min_gen,
        plot_max_gen=plot_max_gen,
        ref_gen=ref_gen,
        g_max_avail=g_max_avail,
    )

    # ---------- phase_diagram, hdist
    phase_diagram = pd_data[g_ref]
    hdist = hdist_data[g_ref]

    # ---------- filtering generations
    if plot_min_gen is not None or plot_max_gen is not None:
        filtered_rslt = rslt_data[(rslt_data['Gen'] >= g_min) & (rslt_data['Gen'] <= g_max)]
        filtered_ids = filtered_rslt.index.values
    else:
        filtered_ids = None

    # ---------- plot
    fig, ax = draw_convex_hull_binary(
        phase_diagram=phase_diagram,
        hdist=hdist,
        filtered_ids=filtered_ids,
        ymax=ymax,
        label_stable=label_stable,
        vmax=vmax,
        bottom_margin=bottom_margin,
        markersize=markersize,
        axis_order=axis_order,
    )

    # ---------- return
    return fig, ax


def plot_convex_hull_ternary(
    plot_min_gen=None,
    plot_max_gen=None,
    ref_gen=None,
    show_max=None,
    label_stable=True,
    vmax=0.2,
    markersize=10,
    axis_order='tlr',
):
    # ---------- load data
    pd_data = pkl_data.load_pd_data()
    hdist_data = pkl_data.load_hdist_data()
    rslt_data = pkl_data.load_rslt()
    rin = pkl_data.load_input()
    atype = rin.atype

    # ---------- generation range
    g_max_avail = max(pd_data.keys())
    g_min, g_max, g_ref = get_generation_range(
        plot_min_gen=plot_min_gen,
        plot_max_gen=plot_max_gen,
        ref_gen=ref_gen,
        g_max_avail=g_max_avail,
    )

    # ---------- phase_diagram, hdist
    phase_diagram = pd_data[g_ref]
    hdist = hdist_data[g_ref]

    # ---------- filtering generations
    if plot_min_gen is not None or plot_max_gen is not None:
        filtered_rslt = rslt_data[(rslt_data['Gen'] >= g_min) & (rslt_data['Gen'] <= g_max)]
        filtered_ids = filtered_rslt.index.values
    else:
        filtered_ids = None

    # ---------- ordering
    ordering = build_ordering(atype, axis_order)

    # ---------- plot
    fig, ax = draw_convex_hull_ternary(
        phase_diagram=phase_diagram,
        hdist=hdist,
        filtered_ids=filtered_ids,
        show_max=show_max,
        label_stable=label_stable,
        vmax=vmax,
        markersize=markersize,
        ordering=ordering,
    )

    # ---------- return
    return fig, ax