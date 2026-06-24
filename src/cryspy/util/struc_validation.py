from pymatgen.core import Structure

from .charge_neutral import is_charge_neutral_nat
from .struc_util import get_nat


def validate_loaded_structures(rin, init_struc_data):
    """
    Validate manually loaded initial structures.

    Parameters
    ----------
    rin : ReadInput
    init_struc_data : dict
        {cid: structure}
    """

    # ---------- data type
    if not isinstance(init_struc_data, dict):
        raise ValueError('init_struc_data must be dict')

    # ---------- number of loaded structures
    n_loaded = len(init_struc_data)
    if rin.algo in ['EA', 'EA-vc']:
        if n_loaded > rin.n_pop:
            raise ValueError(
                f'Number of loaded structures must be less than or equal to n_pop: '
                f'n_pop = {rin.n_pop}, '
                f'actual number of structures = {n_loaded}'
            )
    else:
        if n_loaded > rin.tot_struc:
            raise ValueError(
                f'Number of loaded structures must be less than or equal to tot_struc: '
                f'tot_struc = {rin.tot_struc}, '
                f'actual number of structures = {n_loaded}'
            )

    # ---------- structure IDs
    expected_ids = set(range(n_loaded))
    actual_ids = set(init_struc_data)
    if actual_ids != expected_ids:
        missing_ids = sorted(expected_ids - actual_ids, key=str)
        unexpected_ids = sorted(actual_ids - expected_ids, key=str)
        raise ValueError(
            f'Invalid structure IDs: '
            f'missing IDs = {missing_ids}, '
            f'unexpected IDs = {unexpected_ids}'
        )

    # ---------- structure data
    for cid in range(n_loaded):
        struc = init_struc_data[cid]

        # ------ structure type
        if not isinstance(struc, Structure):
            raise ValueError(
                f'Structure ID {cid}: '
                f'structure data must be pymatgen Structure'
            )

        # ------ composition
        if rin.algo == 'EA-vc':
            _validate_variable_composition_structure(rin, cid, struc)
        else:
            _validate_fixed_composition_structure(rin, cid, struc)


def _validate_variable_composition_structure(rin, cid, struc):
    # ---------- atom types
    allowed_symbols = set(rin.atype)
    actual_symbols = set(struc.symbol_set)
    unexpected_symbols = actual_symbols - allowed_symbols
    if unexpected_symbols:
        raise ValueError(
            f'Structure ID {cid}: '
            f'allowed atom types = {sorted(allowed_symbols)}, '
            f'unexpected atom types = {sorted(unexpected_symbols)}'
        )

    # ---------- number of atoms
    actual_nat = get_nat(struc, rin.atype)
    if sum(actual_nat) == 0:
        raise ValueError(
            f'Structure ID {cid}: structure must contain at least one atom'
        )
    if not all(
            ll <= n_i <= ul
            for n_i, ll, ul in zip(actual_nat, rin.ll_nat, rin.ul_nat)
        ):
        raise ValueError(
            f'Structure ID {cid}: '
            f'expected nat range = [{rin.ll_nat}, {rin.ul_nat}], '
            f'actual nat = {actual_nat}'
        )

    # ---------- composition
    if rin.min_comp is not None or rin.max_comp is not None:
        ntot = sum(actual_nat)
        actual_comp = tuple(n_i / ntot for n_i in actual_nat)
        tol = 1.0e-12
        if (
                rin.min_comp is not None
                and any(
                    comp < lower - tol
                    for comp, lower in zip(actual_comp, rin.min_comp)
                )
            ):
            raise ValueError(
                f'Structure ID {cid}: '
                f'composition is below min_comp: '
                f'min_comp = {rin.min_comp}, '
                f'actual composition = {actual_comp}'
            )
        if (
                rin.max_comp is not None
                and any(
                    comp > upper + tol
                    for comp, upper in zip(actual_comp, rin.max_comp)
                )
            ):
            raise ValueError(
                f'Structure ID {cid}: '
                f'composition is above max_comp: '
                f'max_comp = {rin.max_comp}, '
                f'actual composition = {actual_comp}'
            )

    # ---------- charge neutrality
    if (
            rin.charge is not None
            and not is_charge_neutral_nat(
                actual_nat,
                rin.charge,
                max_cn_grid_points=rin.max_cn_grid_points,
            )
        ):
        raise ValueError(
            f'Structure ID {cid}: '
            f'nat is not charge neutral: actual nat = {actual_nat}'
        )


def _validate_fixed_composition_structure(rin, cid, struc):
    # ---------- atom types
    expected_symbols = {
        atype
        for atype, n_i in zip(rin.atype, rin.nat)
        if n_i > 0
    }
    actual_symbols = set(struc.symbol_set)
    if actual_symbols != expected_symbols:
        raise ValueError(
            f'Structure ID {cid}: '
            f'expected atom types = {sorted(expected_symbols)}, '
            f'actual atom types = {sorted(actual_symbols)}'
        )

    # ---------- number of atoms
    actual_nat = get_nat(struc, rin.atype)
    if actual_nat != rin.nat:
        raise ValueError(
            f'Structure ID {cid}: '
            f'expected nat = {rin.nat}, actual nat = {actual_nat}'
        )