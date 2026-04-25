from ..util.visual_util import (
    draw_convex_hull_binary,
    draw_convex_hull_ternary,
)
from logging import getLogger

import numpy as np
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram


logger = getLogger('cryspy')


def calc_convex_hull(
        atype,
        gen,
        end_point,
        rslt_data,
        nat_data,
        ymax,
        show_max,
        label_stable,
        vmax,
        bottom_margin,
        markersize,
        fig_format,
        emax_ea=None,
        emin_ea=None,
        mpl_draw=True,
        axis_order=None,
        min_comp=None,
        max_comp=None,
        show_comp_window=True,
    ):
    '''
    Input:
        atype (tuple): atom type, e.g. ('Na', 'Cl')
        gen (int): current generation
        end_point (tuple): end points e.g. (0.0  0.0)
        rslt_data (DataFrame): result data
        nat_data (dict): number of atoms of all structures, {ID: (nat1, nat2, ...), ...}
        ymax (float): Binary only. max value of y-axis
        show_max (float): Ternary only. If not None, plots structures with hull distance ≤ show_max
        label_stable (bool): whether to show stable compositions
        vmax (float): max value of colorbar for hull distance
        bottom_margin (float): bottom margin of y-axis
        markersize (int): size of markers
        fig_format (str): format of figure, 'svg', 'png', or 'pdf'.
        emax_ea (float): maximum energy for cutoff
        emin_ea (float): minimum energy for cutoff
        mpl_draw (bool): whether to draw convex hull figure
        axis_order (str): order of axis for binary and ternary phase diagrams
        min_comp (tuple): Minimum composition fractions
        max_comp (tuple): Maximum composition fractions
        show_comp_window (bool): Whether to overlay feasible composition range

    Return:
        hdist [dict]: hull distance of all structures, {ID: distance, ...}
    '''

    # ---------- rslt --> entries
    e_all = rslt_data[rslt_data['Gen'] <= gen]['E_eV_atom'].to_dict()    # energy data under gen
    entries = {}
    for cid, e in e_all.items():
        # ------ np.nan
        if np.isnan(e):
            continue
        # ------ emax_ea
        if emax_ea is not None:
            if e > emax_ea:
                print(f'Eliminate ID {cid} from convex hull: {e} > emax_ea')
                continue
        # ------ emin_ea
        if emin_ea is not None:
            if e < emin_ea:
                print(f'Eliminate ID {cid} from convex hull: {e} < emin_ea')
                continue
        # ------ entry
        composition = "".join(f"{element}{nat_i}" for element, nat_i in zip(atype, nat_data[cid]))
        entries[cid] = ComputedEntry(composition, e*sum(nat_data[cid]), entry_id=cid)

    # ---------- end points
    end_entry_values = [ComputedEntry(element, end_e) for element, end_e in zip(atype, end_point)]

    # ---------- PhaseDiagram and hull distance
    phase_diagram = PhaseDiagram(end_entry_values + list(entries.values()))
    hdist = {cid: phase_diagram.get_e_above_hull(entries[cid]) for cid in entries}

    # ---------- axis order
    if len(atype) == 2:    # ordering not used for binary system
        if axis_order is None:
            axis_order = 'lr'
    if len(atype) == 3:
        if axis_order is None:
            axis_order ='tlr'

    # ---------- draw convex hull
    if mpl_draw:
        if len(atype) == 2:
            fig, _ = draw_convex_hull_binary(
                phase_diagram=phase_diagram,
                hdist=hdist,
                filtered_ids=None,
                ymax=ymax,
                label_stable=label_stable,
                vmax=vmax,
                bottom_margin=bottom_margin,
                markersize=markersize,
                axis_order=axis_order,
                min_comp=min_comp,
                max_comp=max_comp,
                show_comp_window=show_comp_window,
            )
            fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.{fig_format}')
        elif len(atype) == 3:
            fig, _ = draw_convex_hull_ternary(
                atype=atype,
                phase_diagram=phase_diagram,
                hdist=hdist,
                filtered_ids=None,
                show_max=show_max,
                label_stable=label_stable,
                vmax=vmax,
                markersize=markersize,
                axis_order=axis_order,
                min_comp=min_comp,
                max_comp=max_comp,
                show_comp_window=show_comp_window,
            )
            fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.{fig_format}')

    # ---------- return
    return phase_diagram, hdist
