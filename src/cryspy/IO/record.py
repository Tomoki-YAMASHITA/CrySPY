from logging import getLogger

from .out_results import out_rslt
from .pkl_data import save_opt_struc, save_rslt
from ..util.struc_util import out_poscar, out_cif


logger = getLogger('cryspy')


def record_init():
    pass


def record_opt(
        rin,
        cid,
        init_struc_data,
        opt_struc_data,
        rslt_data,
        opt_struc,
        energy,
        magmom,
        check_opt,
        ef=None,
        nat=None,
        n_selection=None,
        gen=None
    ):
    '''
    save opt_sturc_data, rslt_data, and out opt_xxx, rslt_xxx
    '''
    # ---------- get initial spg info
    try:
        spg_sym, spg_num = init_struc_data[
            cid].get_space_group_info(symprec=rin.symprec)
    except TypeError:
        spg_num = 0
        spg_sym = None

    # ---------- success
    if opt_struc is not None:
        # ------ get opt spg info
        try:
            spg_sym_opt, spg_num_opt = opt_struc.get_space_group_info(
                symprec=rin.symprec)
        except:
            spg_num_opt = 0
            spg_sym_opt = None
        # ------ out opt_struc
        out_poscar({cid:opt_struc}, './data/opt_POSCARS')
        try:
            out_cif(opt_struc, cid,
                    './data/opt_CIFS.cif', rin.symprec)
        except TypeError:
            logger.warning('failed to write opt_CIF')
    # ---------- error
    else:
        spg_num_opt = 0
        spg_sym_opt = None

    # ---------- register opt_struc
    opt_struc_data[cid] = opt_struc
    save_opt_struc(opt_struc_data)

    # ---------- register rslt
    if rin.algo in ['RS', 'LAQA']:
        rslt_data.loc[cid] = [spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo == 'BO':
        rslt_data.loc[cid] = [n_selection,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo in ['EA']:
        rslt_data.loc[cid] = [gen,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, magmom, check_opt]
    elif rin.algo in ['EA-vc']:
        rslt_data.loc[cid] = [gen,
                              spg_num, spg_sym, spg_num_opt, spg_sym_opt,
                              energy, ef, nat, magmom, check_opt]
    save_rslt(rslt_data)
    if rin.algo != 'EA-vc':
        out_rslt(rslt_data)
    else:    # EA-vc
        out_rslt(rslt_data, order_ef=True)

    # ---------- return
    return opt_struc_data, rslt_data