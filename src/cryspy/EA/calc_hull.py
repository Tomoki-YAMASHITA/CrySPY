from ..util.visual_util import draw_convex_hull_binary, draw_convex_hull_ternary

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
        show_max,
        label_stable,
        vmax,
        bottom_margin,
        fig_format,
        emax_ea=None,
        emin_ea=None,
        mpl_draw=True,
    ):
    '''
    Input:
        atype (tuple): atom type, e.g. ('Na', 'Cl')
        gen (int): current generation
        end_point (tuple): end points e.g. (0.0  0.0)
        rslt_data (DataFrame): result data
        nat_data (dict): number of atoms of all structures, {ID: (nat1, nat2, ...), ...}
        show_max (float): max value of y-axis (binary) or hull distance (ternary)
        label_stable (bool): whether to show stable compositions
        vmax (float): max value of colorbar for hull distance
        bottom_margin (float): bottom margin of y-axis
        fig_format (str): format of figure, 'svg', 'png', or 'pdf'.
        emax_ea (float): maximum energy for cutoff
        emin_ea (float): minimum energy for cutoff

    Return:
        hdist [dict]: hull distance of all structures, {ID: distance, ...}
    '''

    # ---------- current generation
    c_rslt = rslt_data[rslt_data['Gen'] == gen]
    cgen_ids = c_rslt.index.values    # current IDs [array]

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
    pd = PhaseDiagram(end_entry_values + list(entries.values()))
    hdist = {cid: pd.get_e_above_hull(entries[cid]) for cid in entries}

    # ---------- draw convex hull
    if mpl_draw:
        if len(atype) == 2:
            fig, _ = draw_convex_hull_binary(pd, hdist, cgen_ids, show_max, label_stable, vmax, bottom_margin)
            fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.{fig_format}', bbox_inches='tight')
        elif len(atype) == 3:
            fig, _ = draw_convex_hull_ternary(pd, hdist, cgen_ids, show_max, label_stable, vmax)
            fig.savefig(f'./data/convex_hull/conv_hull_gen_{gen}.{fig_format}', bbox_inches='tight')

    # ---------- return
    return pd, hdist
