#!/usr/bin/env python3
'''
Energy plot script for high-throughput mode
'''

import argparse
import configparser
from contextlib import closing
from pathlib import Path
import sqlite3

from cryspy.IO.read_input import ReadInput
from cryspy.high_throughput.db.ea_db import (
    select_ea_hdist,
    select_ea_composition_bounds,
    select_ea_phase_diagram,
    select_latest_ea_phase_diagram_generation,
)
from cryspy.high_throughput.db.record_db import (
    select_energy_results,
    select_ids_in_generation_range,
)
from cryspy.high_throughput.db.sqlite_db import connect_db
from cryspy.util.visual_util import (
    draw_convex_hull_binary,
    draw_convex_hull_ternary,
    get_compound_comp_window,
    get_generation_range,
    plot_energy_EA,
    plot_energy_RS,
    prepare_compound_comp_window,
    prepare_compound_hull_data,
    prepare_EA_data,
    save_composition_window,
    save_compound_composition_window,
)


DEFAULT_DB_PATH = 'data/db_data/rslt_data.db'


def load_energy_data(conn: sqlite3.Connection):
    """Load completed results for the shared energy plot functions."""
    import pandas as pd

    # ---------- load results
    rslt_data = pd.DataFrame(
        select_energy_results(conn),
        columns=['ID', 'Gen', 'E_eV_atom'],
    ).set_index('ID')
    rslt_data['E_eV_atom'] = pd.to_numeric(rslt_data['E_eV_atom'])
    rslt_data['E_eV_atom'] = rslt_data['E_eV_atom'].replace(
        [float('inf'), -float('inf')],
        float('nan'),
    )

    # ---------- return
    return rslt_data


def load_hull_data(conn: sqlite3.Connection, rin):
    """Load a saved convex hull and select generations for plotting."""
    # ---------- latest saved generation
    latest_generation = select_latest_ea_phase_diagram_generation(conn)
    if latest_generation is None:
        raise ValueError('No saved phase diagram was found')

    # ---------- generation range
    g_min, g_max, g_ref = get_generation_range(
        ref_gen=rin.ref_gen,
        plot_min_gen=rin.plot_min_gen,
        plot_max_gen=rin.plot_max_gen,
        g_max_avail=latest_generation,
    )
    if g_min > g_max:
        raise ValueError('No generations remain in the selected plot range')

    # ---------- phase diagram
    phase_diagram = select_ea_phase_diagram(conn, g_ref)
    if phase_diagram is None:
        raise ValueError(
            f'No saved phase diagram was found for generation {g_ref}'
        )

    # ---------- hull distances
    hdist = select_ea_hdist(conn, g_ref)
    missing_ids = {
        entry.entry_id for entry in phase_diagram.all_entries
        if entry.entry_id is not None
    } - hdist.keys()
    if missing_ids:
        raise ValueError(
            f'Hull distances are missing for generation {g_ref}: '
            f'IDs {sorted(missing_ids)}'
        )

    # ---------- filtering generations
    if rin.plot_min_gen is not None or rin.plot_max_gen is not None:
        filtered_ids = select_ids_in_generation_range(
            conn,
            g_min,
            g_max,
        )
    else:
        filtered_ids = None

    # ---------- composition window
    g_ref_comp = rin.ref_gen_comp if rin.ref_gen_comp is not None else g_ref
    if g_ref_comp > latest_generation:
        raise ValueError(
            f'ref_gen_comp = {g_ref_comp} is larger than the latest generation '
            f'(latest = {latest_generation})'
        )
    comp_bounds = select_ea_composition_bounds(conn, g_ref_comp)
    if comp_bounds is None:
        raise ValueError(f'Generation {g_ref_comp} was not found in ea_info')
    min_comp, max_comp = comp_bounds

    # ---------- return
    return (
        phase_diagram, hdist, min_comp, max_comp,
        filtered_ids, g_ref, g_min, g_max,
    )


def plot_RS(rin, rslt_data) -> Path:
    """Plot RS energies using the shared plotting function."""
    # ---------- check results
    if rslt_data.empty or not rslt_data['E_eV_atom'].notna().any():
        raise ValueError('No finite energies were found for plotting')

    # ---------- plot
    fig, ax = plot_energy_RS(rslt_data, rin.ymax, rin.markersize)

    # ---------- save figure
    output_path = Path('data/energy_plot') / (
        f'energy_RS_{len(rslt_data)}.{rin.fig_format}'
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    return output_path


def plot_EA(rin, rslt_data) -> Path:
    """Plot EA energies with the normal-mode generation selection."""
    # ---------- check results
    if rslt_data.empty or rslt_data['Gen'].isna().any():
        raise ValueError('No valid EA generation results were found')

    # ---------- prepare data
    filtered_rslt, max_indx, g_ref, g_min, g_max = prepare_EA_data(
        rslt_data=rslt_data,
        ref_gen=rin.ref_gen,
        plot_min_gen=rin.plot_min_gen,
        plot_max_gen=rin.plot_max_gen,
    )
    if filtered_rslt.empty or not filtered_rslt['E_eV_atom'].notna().any():
        raise ValueError('No finite energies were found in the selected generations')

    # ---------- plot
    fig, ax = plot_energy_EA(
        filtered_rslt,
        rin.ymax,
        rin.markersize,
        max_indx=max_indx,
    )

    # ---------- save figure
    if g_min == 1 and g_max == g_ref:
        filename = f'energy_EA_gen_{g_ref}.{rin.fig_format}'
    else:
        filename = (
            f'energy_EA_gen_{g_min}-{g_max}'
            f'_ref_{g_ref}.{rin.fig_format}'
        )
    output_path = Path('data/energy_plot') / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    return output_path


def plot_EA_vc(rin, hull_data) -> Path:
    """Plot a saved binary or ternary convex hull with optional composition windows."""
    # ---------- prepare data
    (
        phase_diagram,
        hdist,
        min_comp,
        max_comp,
        filtered_ids,
        g_ref,
        g_min,
        g_max,
    ) = hull_data
    plot_components = rin.atype
    comp_window = None
    if rin.terminals is not None:
        phase_diagram, hdist, plot_components = prepare_compound_hull_data(
            phase_diagram,
            rin.terminals,
        )
        if rin.show_comp_window:
            comp_window = get_compound_comp_window(
                atype=rin.atype,
                compound_phase_diagram=phase_diagram,
                min_comp=min_comp,
                max_comp=max_comp,
            )
        min_comp = None
        max_comp = None
    elif not rin.show_comp_window:
        min_comp = None
        max_comp = None

    # ---------- axis order
    axis_order = (
        rin.terminal_axis_order
        if rin.terminals is not None
        else rin.axis_order
    )
    if axis_order is None:
        axis_order = 'lr' if len(plot_components) == 2 else 'tlr'

    # ---------- plot
    if len(plot_components) == 2:
        fig, ax = draw_convex_hull_binary(
            phase_diagram=phase_diagram,
            hdist=hdist,
            filtered_ids=filtered_ids,
            ymax=rin.ymax,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            bottom_margin=rin.bottom_margin,
            markersize=rin.markersize,
            axis_order=axis_order,
            min_comp=min_comp,
            max_comp=max_comp,
            comp_window=comp_window,
        )
    elif len(plot_components) == 3:
        fig, ax = draw_convex_hull_ternary(
            atype=plot_components,
            phase_diagram=phase_diagram,
            hdist=hdist,
            filtered_ids=filtered_ids,
            show_max=rin.show_max,
            label_stable=rin.label_stable,
            vmax=rin.vmax,
            markersize=rin.markersize,
            axis_order=axis_order,
            min_comp=min_comp,
            max_comp=max_comp,
            comp_window=comp_window,
        )
    else:
        raise ValueError('Convex hull plots require two or three plot components')

    # ---------- save figure
    terminal_suffix = (
        f'_terminals_{"-".join(plot_components)}'
        if rin.terminals is not None
        else ''
    )
    if g_min == 1 and g_max == g_ref:
        filename = f'conv_hull_gen_{g_ref}{terminal_suffix}.{rin.fig_format}'
    else:
        filename = (
            f'conv_hull_gen_{g_min}-{g_max}'
            f'_ref_{g_ref}{terminal_suffix}.{rin.fig_format}'
        )
    output_path = Path('data/convex_hull') / filename
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    return output_path


def plot_composition_window_input(rin) -> Path:
    """Plot the composition window from cryspy.in."""

    # ---------- check input
    if rin.algo != 'EA-vc':
        raise ValueError('Composition window plot is available only for EA-vc.')
    if rin.min_comp is None and rin.max_comp is None:
        raise ValueError('min_comp and max_comp are not set in cryspy.in.')

    # ---------- elemental composition window
    if rin.terminals is None:
        if len(rin.atype) not in (2, 3):
            raise ValueError(
                f'len(rin.atype) = {len(rin.atype)} is not supported '
                'for composition window plot'
            )
        save_composition_window(
            atype=rin.atype,
            gen='input',
            min_comp=rin.min_comp,
            max_comp=rin.max_comp,
            fig_format=rin.fig_format,
            axis_order=rin.axis_order,
        )
        return Path(f'data/convex_hull/composition_window_input.{rin.fig_format}')

    # ---------- compound composition window
    terminal_labels, comp_window = prepare_compound_comp_window(
        atype=rin.atype,
        terminals=rin.terminals,
        min_comp=rin.min_comp,
        max_comp=rin.max_comp,
    )
    return Path(save_compound_composition_window(
        terminal_labels=terminal_labels,
        gen='input',
        comp_window=comp_window,
        fig_format=rin.fig_format,
        axis_order=rin.terminal_axis_order,
    ))


def main():
    """Plot CrySPY high-throughput results from SQLite."""

    # ---------- argparse
    parser = argparse.ArgumentParser(
        description='Plot CrySPY high-throughput energies and convex hulls',
        epilog='Plot settings are read from [visual] in cryspy.in.',
    )
    parser.add_argument(
        'dbfile',
        nargs='?',
        default=DEFAULT_DB_PATH,
        help=f'SQLite database (default: {DEFAULT_DB_PATH})',
    )
    parser.add_argument(
        '-c',
        '--composition-window',
        action='store_true',
        help='plot composition window from cryspy.in',
    )
    args = parser.parse_args()

    try:
        # ---------- read input
        rin = ReadInput(ht=True)
        if rin.algo not in ('RS', 'EA', 'EA-vc'):
            raise ValueError('cryspy-ht-Eplot supports RS, EA, and EA-vc')

        # ---------- composition window
        if args.composition_window:
            output_path = plot_composition_window_input(rin)
            print(f'Output: {output_path}')
            return

        # ---------- check database
        db_path = Path(args.dbfile)
        if not db_path.is_file():
            raise ValueError(f'{db_path} does not exist')

        # ---------- database snapshot
        with closing(connect_db(str(db_path))) as conn:
            conn.execute('BEGIN')
            if rin.algo == 'EA-vc':
                hull_data = load_hull_data(conn, rin)
            else:
                rslt_data = load_energy_data(conn)

        # ---------- plot
        if rin.algo == 'RS':
            output_path = plot_RS(rin, rslt_data)
        elif rin.algo == 'EA':
            output_path = plot_EA(rin, rslt_data)
        else:
            output_path = plot_EA_vc(rin, hull_data)

        # ---------- output
        print(f'Database: {db_path}')
        print(f'Output: {output_path}')

    except (OSError, ValueError, sqlite3.Error, configparser.Error) as e:
        parser.error(str(e))


if __name__ == '__main__':
    # ---------- main
    main()
