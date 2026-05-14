from logging import getLogger

import numpy as np

from .constants import MAX_CN_GRID_POINTS

# ---------- import later
#from ..IO.pkl_data import load_cn_comb_data, save_cn_comb_data
#from .struc_util import get_feasible_composition, precompute_feasible_N
#from .struc_util import sample_nat_from_feasible_N


logger = getLogger('cryspy')


def prepare_cn_data(
        ll_nat,
        ul_nat,
        charge,
        cn_mode='auto',
        max_cn_grid_points=MAX_CN_GRID_POINTS,
        use_pkl=True,
    ):
    if cn_mode not in ['auto', 'enumerate', 'sample']:
        raise ValueError('cn_mode must be auto, enumerate, or sample')

    # ---------- use pkl data if it exists and matches the current parameters
    if use_pkl:
        from ..IO.pkl_data import load_cn_comb_data
        try:
            cn_data = load_cn_comb_data()
        except FileNotFoundError:
            cn_data = None

        if cn_data is not None:
            if not isinstance(cn_data, dict):
                raise ValueError(
                    'cn_comb_data.pkl has an old unsupported format. '
                    'The data format was changed in CrySPY 1.5.0. '
                    'Please use a consistent CrySPY version.'
                )

            if (
                    cn_data.get('ll_nat') == ll_nat
                    and cn_data.get('ul_nat') == ul_nat
                    and cn_data.get('charge') == charge
                    and cn_data.get('cn_mode') == cn_mode
                    and cn_data.get('max_cn_grid_points') == max_cn_grid_points
                ):

                return cn_data

    if cn_mode == 'sample':
        logger.info('# ---------- Prepare charge-neutral data')
        logger.info('Charge-neutral mode: sample')
        cn_data = {
            'mode': 'sample',
            'cn_mode': cn_mode,
            'll_nat': ll_nat,
            'ul_nat': ul_nat,
            'charge': charge,
            'max_cn_grid_points': max_cn_grid_points,
            'cn_comb': None,
        }
    else:
        try:
            logger.info('# ---------- Prepare charge-neutral data')
            logger.info('Charge-neutral mode: enumerate')
            cn_comb = calc_cn_comb(
                ll_nat,
                ul_nat,
                charge,
                max_cn_grid_points=max_cn_grid_points,
            )
            cn_data = {
                'mode': 'enumerate',
                'cn_mode': cn_mode,
                'll_nat': ll_nat,
                'ul_nat': ul_nat,
                'charge': charge,
                'max_cn_grid_points': max_cn_grid_points,
                'cn_comb': cn_comb,
            }
        except ValueError as e:
            if 'Too many charge-neutral enumeration grid points' not in str(e):
                raise
            if cn_mode == 'enumerate':
                raise
            logger.info(str(e))
            logger.info('Charge-neutral mode: sample')
            cn_data = {
                'mode': 'sample',
                'cn_mode': cn_mode,
                'll_nat': ll_nat,
                'ul_nat': ul_nat,
                'charge': charge,
                'max_cn_grid_points': max_cn_grid_points,
                'cn_comb': None,
            }

    # ---------- save data to pkl
    if use_pkl:
        from ..IO.pkl_data import save_cn_comb_data
        save_cn_comb_data(cn_data)
        if cn_data['mode'] == 'enumerate':
            logger.info(f'Charge neutral combinations saved: {len(cn_data["cn_comb"])}')
        else:
            logger.info('Charge neutral sampler data saved')

    # ---------- return
    return cn_data


def calc_cn_comb(ll_nat, ul_nat, charge, max_cn_grid_points=MAX_CN_GRID_POINTS):
    """
    Calculate charge-neutral combinations of atoms

    # ---------- args
    ll_nat (tuple): lower limit of the number of atoms
    ul_nat (tuple): upper limit of the number of atoms
    charge (tuple): charge of each atom type
    max_cn_grid_points (int): maximum number of independent atom-count grid points
                              for exhaustive enumeration

    # ---------- retrun
    cn_comb (np.ndarray): array of charge-neutral combinations

    # ---------- example
    e.g. Na-Cl
    ll_nat = (0, 0)
    ul_nat = (4, 4)
    charge = (1, -1)

    cn_comb = array([[1, 1],
                     [2, 2],
                     [3, 3],
                     [4, 4]])

    # ---------- multi-valence example
    e.g. Fe-O
    ll_nat = (4, 4)
    ul_nat = (10, 10)
    charge = ((2, 3), -2)   # (Fe2+, Fe3+), O2-

    Expanded representation:
    ll_ext = (0, 0, 4)      # multivalent: no ll per valence; ll applies to total only
    ul_ext = (10, 10, 10)
    charge_ext = (2, 3, -2)
    """

    # ---------- convert bounds to np.array
    ll = np.array(ll_nat, dtype=int)
    ul = np.array(ul_nat, dtype=int)
    k_orig = len(charge)

    # ---------- case 1: all charges are single-valence (int)
    if all(isinstance(c, int) for c in charge):
        charges = np.array(charge, dtype=int)
        cn_comb = _calc_cn_comb_single(
            ll,
            ul,
            charges,
            max_cn_grid_points=max_cn_grid_points,
        )

    # ---------- case 2: some species have multiple valences
    else:
        ll_ext, ul_ext, charges_ext, orig_index = _expand_charge_valence(
            ll,
            ul,
            charge,
        )

        # ------ charge-neutral combinations in the extended space
        cn_ext = _calc_cn_comb_single(
            ll_ext,
            ul_ext,
            charges_ext,
            max_cn_grid_points=max_cn_grid_points,
        )
        if cn_ext.size == 0:
            return np.empty((0, k_orig), dtype=int)

        # ------ aggregate back to original species
        totals = np.zeros((cn_ext.shape[0], k_orig), dtype=int)    # initialized
        for j, p in enumerate(orig_index):
            totals[:, p] += cn_ext[:, j]

        # ------ apply original bounds ll_nat, ul_nat
        mask_ll = np.all(totals >= ll, axis=1)
        mask_ul = np.all(totals <= ul, axis=1)
        final_mask = mask_ll & mask_ul
        if not np.any(final_mask):
            return np.empty((0, k_orig), dtype=int)
        totals_valid = totals[final_mask]

        # ------ remove duplicate compositions
        cn_comb = np.unique(totals_valid, axis=0)

    # ---------- return charge neutral combinations
    return cn_comb


def _calc_cn_comb_single(ll, ul, charges, max_cn_grid_points=MAX_CN_GRID_POINTS):
    """
    Core routine for charge-neutral combinations with single valence per species.

    Parameters
    ----------
    ll : np.ndarray, shape (k,)
        Lower limits of the number of atoms for each species (int).
    ul : np.ndarray, shape (k,)
        Upper limits of the number of atoms for each species (int).
    charges : np.ndarray, shape (k,)
        Charges of each species (single valence only, all ints).
    max_cn_grid_points : int
        Maximum number of independent atom-count grid points for exhaustive enumeration.

    Returns
    -------
    cn_comb : np.ndarray, shape (N, k)
        Array of charge-neutral combinations.
        Each row is (n_1, ..., n_k).
    """

    # ---------- size
    k = len(charges)

    # ---------- choose one "dependent" species whose charge != 0
    #            Select the dependent species: charge with largest absolute value
    nonzero_q_idx = np.where(charges != 0)[0]
    dep_idx = nonzero_q_idx[np.argmax(np.abs(charges[nonzero_q_idx]))]    # dependent index
    indep_idx = [i for i in range(k) if i != dep_idx]
    q_dep = charges[dep_idx]
    q_indep = charges[indep_idx]

    # ---------- check number of combinations before meshgrid
    n_grid_points = 1
    for i in indep_idx:
        n_grid_points *= int(ul[i] - ll[i] + 1)
    if n_grid_points > max_cn_grid_points:
        raise ValueError(
            f'Too many charge-neutral enumeration grid points: {n_grid_points} '
            f'(max_cn_grid_points = {max_cn_grid_points})'
        )

    # ---------- generate combinations for independent variables only
    ranges = [np.arange(ll[i], ul[i] + 1, dtype=np.int32) for i in indep_idx]
    mesh = np.meshgrid(*ranges, indexing='ij')
    comb_indep = np.stack(mesh, axis=-1).reshape(-1, len(indep_idx))  # (M, k-1)

    # ---------- compute dependent variable from charge neutrality
    # q_dep * n_dep + sum(q_indep * n_indep) = 0
    # → n_dep = - sum(q_indep * n_indep) / q_dep
    sum_indep = comb_indep.astype(np.int64) @ q_indep.astype(np.int64)  # (M,)
    divisible_mask = (-sum_indep % q_dep == 0)    # divisible by q_dep
    n_dep = -sum_indep // q_dep

    # ---------- mask
    in_range_mask = (n_dep >= ll[dep_idx]) & (n_dep <= ul[dep_idx])    # ll <= n_dep <= ul
    total_atoms = comb_indep.sum(axis=1) + n_dep    # \sum n_indep + n_dep
    nonzero_mask = (total_atoms != 0)    # total atoms != 0
    final_mask = divisible_mask & in_range_mask & nonzero_mask

    # ---------- construct final combinations
    if not np.any(final_mask):
        cn_comb = np.empty((0, k), dtype=int)    # no valid combinations --> empty array
    else:
        valid_indep = comb_indep[final_mask]
        valid_dep = n_dep[final_mask]
        # ------ reconstruct full combinations
        cn_comb = np.zeros((len(valid_dep), k), dtype=int)    # vacant array
        cn_comb[:, dep_idx] = valid_dep
        for j, i in enumerate(indep_idx):
            cn_comb[:, i] = valid_indep[:, j]

    # ---------- return
    return cn_comb


def filter_cn_comb_comp(cn_comb, min_comp=None, max_comp=None, tol=1e-12):
    """
    Filter charge-neutral combinations by composition window.
    """
    if cn_comb is None or len(cn_comb) == 0:
        return np.empty((0, 0), dtype=int)

    comp = cn_comb / cn_comb.sum(axis=1, keepdims=True)

    mask = np.ones(len(cn_comb), dtype=bool)    # all True initially
    if min_comp is not None:
        min_comp = np.asarray(min_comp, dtype=float)
        mask &= np.all(comp >= (min_comp - tol), axis=1)
    if max_comp is not None:
        max_comp = np.asarray(max_comp, dtype=float)
        mask &= np.all(comp <= (max_comp + tol), axis=1)

    return cn_comb[mask]


def _expand_charge_valence(ll, ul, charge):
    """
    Expand each valence to a pseudo-species.
    """
    valence_list = []    # charges in extended space
    orig_index = []      # mapping: extended index -> original species index
    ll_ext_list = []     # lower bounds in extended space

    for i, c in enumerate(charge):
        if isinstance(c, int):
            # ---------- single-valence species: lower bound can be kept
            valence_list.append(c)
            orig_index.append(i)
            ll_ext_list.append(ll[i])
        else:
            # ---------- multivalent species: each valence may be absent
            for v in c:
                valence_list.append(v)
                orig_index.append(i)
                ll_ext_list.append(0)   # no ll per valence; ll applies to total only

    # ---------- convert to np.array
    charges_ext = np.array(valence_list, dtype=int)
    ll_ext = np.array(ll_ext_list, dtype=int)
    ul_ext = np.array([ul[i] for i in orig_index], dtype=int)

    # ---------- return
    return ll_ext, ul_ext, charges_ext, orig_index


def sample_cn_nat(
        ll_nat,
        ul_nat,
        charge,
        max_trials=100000,
        rng=None,
    ):
    """
    Sample one charge-neutral nat.
    """
    if rng is None:
        rng = np.random.default_rng()

    # ---------- convert bounds to np.array
    ll = np.array(ll_nat, dtype=int)
    ul = np.array(ul_nat, dtype=int)
    k_orig = len(charge)

    # ---------- case 1: all charges are single-valence (int)
    if all(isinstance(c, int) for c in charge):
        charges = np.array(charge, dtype=int)
        for _ in range(max_trials):
            cn_comb = _sample_cn_comb_single(ll, ul, charges, rng)
            if cn_comb is not None:
                return tuple(int(n) for n in cn_comb[0])

    # ---------- case 2: some species have multiple valences
    else:
        ll_ext, ul_ext, charges_ext, orig_index = _expand_charge_valence(
            ll,
            ul,
            charge,
        )

        for _ in range(max_trials):
            cn_ext = _sample_cn_comb_single(ll_ext, ul_ext, charges_ext, rng)
            if cn_ext is None:
                continue

            # ------ aggregate back to original species
            nat = np.zeros(k_orig, dtype=int)    # initialized
            for j, p in enumerate(orig_index):
                nat[p] += cn_ext[0, j]

            # ------ apply original bounds ll_nat, ul_nat
            if np.any(nat < ll) or np.any(nat > ul):
                continue

            return tuple(int(n) for n in nat)

    raise ValueError(
        f'Failed to sample a charge-neutral nat in {max_trials} trials. '
        'Please check ll_nat, ul_nat, and charge.'
    )


def _sample_cn_comb_single(ll, ul, charges, rng):
    """
    Sample one charge-neutral combination in single-valence/extended space.
    Positive-negative dependent pairs are tried in random order.
    """
    # ---------- size
    k = len(charges)

    # ---------- positive and negative species
    pos_idx = np.where(charges > 0)[0]
    neg_idx = np.where(charges < 0)[0]
    if len(pos_idx) == 0 or len(neg_idx) == 0:
        return None

    # ---------- make positive-negative pairs and randomize order
    pairs = [(p, m) for p in pos_idx for m in neg_idx]
    pair_order = rng.permutation(len(pairs))

    # ---------- try each pair
    for i_pair in pair_order:
        p_idx, m_idx = pairs[i_pair]
        pair_idx = [p_idx, m_idx]
        indep_idx = [i for i in range(k) if i not in pair_idx]
        q_p = charges[p_idx]
        q_m = charges[m_idx]

        # ------ sample independent variables
        cn_comb = np.zeros((1, k), dtype=int)
        if indep_idx:
            cn_comb[0, indep_idx] = rng.integers(ll[indep_idx], ul[indep_idx] + 1)

        # ------ compute target charge for dependent pair
        # q_p * n_p + q_m * n_m + sum(q_i * n_i) = 0
        # --> q_p * n_p + q_m * n_m = -sum(q_i * n_i)
        target = 0
        if indep_idx:
            target = -int(cn_comb[0, indep_idx] @ charges[indep_idx])

        # ------ find bounded integer solutions for the dependent pair
        candidates = []
        for n_p in range(ll[p_idx], ul[p_idx] + 1):
            numerator = target - q_p * n_p
            if numerator % q_m != 0:
                continue
            n_m = numerator // q_m
            if ll[m_idx] <= n_m <= ul[m_idx]:
                candidates.append((n_p, n_m))
        if len(candidates) == 0:
            continue

        # ------ choose one solution
        n_p, n_m = candidates[rng.integers(len(candidates))]
        cn_comb[0, p_idx] = n_p
        cn_comb[0, m_idx] = n_m

        # ------ total atoms must not be zero
        if cn_comb.sum() == 0:
            continue

        # ------ valid combination found
        return cn_comb

    # ---------- if no pair yields a valid combination, return None
    return None


def sample_cn_nat_comp(
        ll_nat,
        ul_nat,
        charge,
        min_comp,
        max_comp,
        feasible_N=None,
        max_trials=100000,
        rng=None,
    ):
    """
    Sample one charge-neutral nat under composition constraints.
    """
    # ---------- rng
    if rng is None:
        rng = np.random.default_rng()

    from .struc_util import get_feasible_composition, precompute_feasible_N

    # ---------- prepare feasible total atom counts
    if feasible_N is None:
        feasible_comp = get_feasible_composition(min_comp, max_comp)
        if feasible_comp is None:
            raise ValueError('No feasible composition exists for min_comp and max_comp.')

        # ------ pre-check charge neutrality in composition space
        charge_range = _charge_range_under_comp(feasible_comp, charge)
        if charge_range is None:
            raise ValueError(
                'No charge-neutral composition exists within min_comp and max_comp.'
            )
        feasible_N = precompute_feasible_N(ll_nat, ul_nat, feasible_comp)

    # ---------- check feasible total atom counts
    if len(feasible_N) == 0:
        raise ValueError('No feasible total atom count exists for composition constraints.')

    # ---------- sample charge-neutral nat under composition constraints
    k_orig = len(charge)
    for _ in range(max_trials):
        # ------ choose a feasible total atom count
        N, lower, upper = feasible_N[rng.integers(len(feasible_N))]

        # ------ single-valence case
        if all(isinstance(c, int) for c in charge):
            charges = np.array(charge, dtype=int)
            cn_comb = _sample_cn_comb_single_N(lower, upper, charges, N, rng)
            if cn_comb is not None:
                return tuple(int(n) for n in cn_comb[0])

        # ------ multi-valence case: expand to extended space, sample, then aggregate back
        else:
            ll_ext, ul_ext, charges_ext, orig_index = _expand_charge_valence(
                lower,
                upper,
                charge,
            )
            cn_ext = _sample_cn_comb_single_N(ll_ext, ul_ext, charges_ext, N, rng)
            if cn_ext is None:
                continue

            # -- aggregate back to original species
            nat = np.zeros(k_orig, dtype=int)    # initialized
            for j, p in enumerate(orig_index):
                nat[p] += cn_ext[0, j]

            # -- apply original bounds for this N
            if np.any(nat < lower) or np.any(nat > upper):
                continue

            return tuple(int(n) for n in nat)

    raise ValueError(
        f'Failed to sample a charge-neutral nat with composition constraints '
        f'in {max_trials} trials.'
    )


def _sample_cn_comb_single_N(ll, ul, charges, N, rng):
    """
    Sample one charge-neutral combination with fixed total atom count N.
    Positive-negative dependent pairs are tried in random order.
    """
    from .struc_util import sample_nat_from_feasible_N

    # ---------- size
    k = len(charges)

    # ---------- positive and negative species
    pos_idx = np.where(charges > 0)[0]
    neg_idx = np.where(charges < 0)[0]
    if len(pos_idx) == 0 or len(neg_idx) == 0:
        return None

    # ---------- make positive-negative pairs and prioritize wider bounds
    width = (ul - ll).astype(int)
    pairs = [(p, m) for p in pos_idx for m in neg_idx]
    tie_break = rng.random(len(pairs))
    pair_order = sorted(
        range(len(pairs)),
        key=lambda i: (-(width[pairs[i][0]] + width[pairs[i][1]]), tie_break[i]),
    )

    # ---------- try each pair
    for i_pair in pair_order:
        p_idx, m_idx = pairs[i_pair]
        pair_idx = [p_idx, m_idx]
        indep_idx = sorted(
            [i for i in range(k) if i not in pair_idx],
            key=lambda i: width[i],
        )
        cn_comb = np.zeros((1, k), dtype=int)

        # ------ sample independent variables with fixed total atom count
        if indep_idx:
            lower_indep = ll[indep_idx]
            upper_indep = ul[indep_idx]

            # -- feasible range of total atoms for independent variables
            N_indep_min = max(
                int(lower_indep.sum()),
                int(N - ul[p_idx] - ul[m_idx]),
            )
            N_indep_max = min(
                int(upper_indep.sum()),
                int(N - ll[p_idx] - ll[m_idx]),
            )
            if N_indep_min > N_indep_max:
                continue

            N_indep = rng.integers(N_indep_min, N_indep_max + 1)
            nat_indep, _ = sample_nat_from_feasible_N(
                [(N_indep, lower_indep, upper_indep)],
                rng,
            )
            cn_comb[0, indep_idx] = nat_indep

        # ------ compute remaining atoms and charge for the dependent pair
        sum_indep = int(cn_comb[0, indep_idx].sum()) if indep_idx else 0
        charge_indep = int(cn_comb[0, indep_idx] @ charges[indep_idx]) if indep_idx else 0
        S_pair = int(N - sum_indep)
        Q_pair = -charge_indep
        q_p = int(charges[p_idx])
        q_m = int(charges[m_idx])

        # ------ solve the two-variable linear equations
        numerator = Q_pair - q_m * S_pair
        denominator = q_p - q_m
        if numerator % denominator != 0:
            continue
        n_p = numerator // denominator
        n_m = S_pair - n_p

        # ------ check bounds for the dependent pair
        if not (ll[p_idx] <= n_p <= ul[p_idx]):
            continue
        if not (ll[m_idx] <= n_m <= ul[m_idx]):
            continue
        cn_comb[0, p_idx] = n_p
        cn_comb[0, m_idx] = n_m

        # ------ final check
        if cn_comb.sum() != N:
            continue
        if int(cn_comb[0] @ charges) != 0:
            continue

        # ------ valid combination found
        return cn_comb

    # ---------- if no pair yields a valid combination, return None
    return None


def _charge_range_under_comp(feasible_comp, charge, tol=1e-12):
    """
    Return min/max charge per atom under composition constraints.
    Multi-valence species are treated by their minimum/maximum valence.
    """
    # ---------- convert composition bounds to np.array
    low = np.array([c[0] for c in feasible_comp], dtype=float)
    high = np.array([c[1] for c in feasible_comp], dtype=float)

    # ---------- valence range for each species
    q_low = np.array(
        [min(c) if isinstance(c, tuple) else c for c in charge],
        dtype=float,
    )
    q_high = np.array(
        [max(c) if isinstance(c, tuple) else c for c in charge],
        dtype=float,
    )

    def optimize(q, order):
        # ------ start from the lower composition bounds
        comp = low.copy()
        remaining = 1.0 - comp.sum()

        # ------ distribute remaining composition greedily
        for i in order:
            if remaining <= tol:
                break
            add = min(remaining, high[i] - low[i])
            comp[i] += add
            remaining -= add

        # ------ return charge per atom
        return float(comp @ q)

    # ---------- minimize and maximize charge per atom
    q_min = optimize(q_low, np.argsort(q_low))
    q_max = optimize(q_high, np.argsort(-q_high))

    # ---------- check whether zero charge is included
    if q_min > tol or q_max < -tol:
        return None

    # ---------- return
    return q_min, q_max


def is_charge_neutral_nat(nat, charge, max_cn_grid_points=MAX_CN_GRID_POINTS):
    """
    Check whether nat can be charge neutral.
    """
    # ---------- single-valence case
    if all(isinstance(c, int) for c in charge):
        charges = np.array(charge, dtype=int)
        return int(np.array(nat, dtype=int) @ charges) == 0

    # ---------- convert to np.array
    nat = np.array(nat, dtype=int)

    # ---------- expand multi-valence species
    ll_ext = []
    ul_ext = []
    charge_ext = []
    orig_index = []

    for i, (n, c) in enumerate(zip(nat, charge)):
        if isinstance(c, int):
            ll_ext.append(int(n))
            ul_ext.append(int(n))
            charge_ext.append(c)
            orig_index.append(i)
        else:
            for v in c:
                ll_ext.append(0)
                ul_ext.append(int(n))
                charge_ext.append(v)
                orig_index.append(i)

    ll_ext = np.array(ll_ext, dtype=int)
    ul_ext = np.array(ul_ext, dtype=int)
    charge_ext = np.array(charge_ext, dtype=int)

    # ---------- enumerate charge-neutral combinations in extended space
    cn_ext = _calc_cn_comb_single(
        ll_ext,
        ul_ext,
        charge_ext,
        max_cn_grid_points=max_cn_grid_points,
    )

    # ---------- check whether any extended combination matches the original nat
    if len(cn_ext) == 0:
        return False

    totals = np.zeros((cn_ext.shape[0], len(nat)), dtype=int)
    for j, p in enumerate(orig_index):
        totals[:, p] += cn_ext[:, j]

    return bool(np.any(np.all(totals == nat, axis=1)))

