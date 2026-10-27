from dataclasses import dataclass
from itertools import product
from logging import getLogger
from typing import Optional

import numpy as np
from pymatgen.core import Structure

from ..ea_offspring import (
    GenerationResult,
    OffspringResult,
    ParentData,
    SelectionContext,
    generate_with_parent_attempts,
    make_task_rng,
)
from ...util.struc_util import check_distance, sort_by_atype, get_nat


logger = getLogger('cryspy')


@dataclass(frozen=True)
class SubstitutionContext:
    """Parameters for substitution generation."""

    atype: tuple[str, ...]
    mindist: tuple[tuple[float, ...], ...]
    subs_comb_map: dict[int, list[list[tuple[str, str]]]]
    maxcnt_ea: int


def gen_child(atype, mindist, parent_A, subs_element_list, maxcnt_ea=50, rng=None):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    subs_element_list (list): list of tuples, e.g. [('Li', 'Ca'), ('Ca', 'Cl')]
    maxcnt_ea (int): maximum number of trial in crossover
    rng (np.random.Generator): random number generator

    # ---------- return
    (if success) child (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    cnt = 0

    # ---------- generate structure by substitution
    while True:
        child = parent_A.copy()    # keep original structure
        # ------ available indices for each atom type
        species_list = [str(site.specie) for site in child]
        available_indices = {elem: [i for i, sp in enumerate(species_list) if sp == elem] 
                            for elem in atype}
        # ------ select indices for substitution without duplicates
        subs_indices = []
        used_indices = set()
        for from_elem, to_elem in subs_element_list:
            available = [idx for idx in available_indices[from_elem] if idx not in used_indices]
            if not available:
                logger.error(f'Substitution invariant violated: no available indices for {from_elem}')
                raise RuntimeError(f'Substitution: no available indices for {from_elem}')
            selected_idx = rng.choice(available)
            subs_indices.append(selected_idx)
            used_indices.add(selected_idx)
        # ------ substitution
        for idx, (from_elem, to_elem) in zip(subs_indices, subs_element_list):
            child.replace(idx, to_elem)
        # ------ check mindist
        success, mindist_ij, dist = check_distance(child, atype, mindist)
        if success:
            child = sort_by_atype(child, atype)
            return child
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in substitution: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Substitution: cnt >= maxcnt_ea. Change parent.')
                return None    # change parent


class SubstitutionOffspringGenerator:
    """Generate one offspring by substitution."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: SubstitutionContext,
        max_parent_attempts: int,
    ) -> None:
        # ---------- check input
        if max_parent_attempts < 1:
            raise ValueError(
                'max_parent_attempts must be greater than zero'
            )

        # ---------- context
        self.parent_data = parent_data
        self.selection_context = selection_context
        self.context = context
        self.max_parent_attempts = max_parent_attempts

    def generate(
        self,
        rng: np.random.Generator,
    ) -> GenerationResult:
        """Generate one offspring."""

        # ---------- generate offspring
        return generate_with_parent_attempts(
            operation='substitution',
            selection_context=self.selection_context,
            n_parent=1,
            max_parent_attempts=self.max_parent_attempts,
            generate_attempt=self._generate_attempt,
            rng=rng,
        )

    def _generate_attempt(
        self,
        parent_ids: tuple[int, ...],
        rng: np.random.Generator,
    ) -> Optional[Structure]:
        """Attempt to generate one offspring."""

        # ---------- parent
        parent_id = parent_ids[0]
        parent = self.parent_data.structures[parent_id]

        # ---------- substitution combinations
        subs_comb = self.context.subs_comb_map[parent_id]
        subs_element_list = subs_comb[int(rng.integers(len(subs_comb)))]

        # ---------- generate offspring
        return gen_child(
            atype=self.context.atype,
            mindist=self.context.mindist,
            parent_A=parent,
            subs_element_list=subs_element_list,
            maxcnt_ea=self.context.maxcnt_ea,
            rng=rng,
        )


def gen_substitution_batch(
        generator,
        n_subs,
        id_start,
        symprec,
        atype,
        base_seed,
        rng,
    ):
    """Generate substitution offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_subs):
        task_rng = make_task_rng(base_seed, cid, rng)
        result = generator.generate(task_rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring_data[cid] = result.structure
        parents[cid] = result.parent_ids
        operation[cid] = 'substitution'
        try:
            spg_sym, spg_num = result.structure.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(result.structure, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} by substitution.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation


def get_subs_comb(atype, ll_nat, ul_nat, subs_max, parent_nat, cn_set=None, charge=None):
    """
    Generate candidate substitution patterns.

    This function returns all feasible substitution patterns under the
    lower-bound, upper-bound, and maximum-substitution constraints.

    Parameters
    ----------
    atype : tuple[str]
        Atom types.
    ll_nat : tuple[int]
        Per-element lower bounds of atom counts.
    ul_nat : tuple[int]
        Per-element upper bounds of atom counts.
    subs_max : int
        Maximum total number of atoms to substitute in one operation.
    parent_nat : tuple[int]
        nat of the parent structure.
    cn_set : set[tuple[int, ...]]
        Set of charge-neutral nat tuples from precomputed cn_comb.
    charge : tuple
        Charge of each atom type.

    Returns
    -------
    subs_comb : list[list[tuple[str, str]]]
        List of substitution patterns.
        Each pattern is a list of (from_elem, to_elem) pairs.

    Notes
    -----
    If cn_set is given, charge neutrality is checked by membership in cn_set.
    If cn_set is None and charge is given, charge neutrality is checked from
    new_nat after substitution.
    """
    # ---------- charge-neutral checker
    if charge is not None and cn_set is None:
        from ...util.charge_neutral import is_charge_neutral_nat

    # ---------- initialize
    subs_comb = []

    # ---------- maximum number of subs for each element pair (no constraints)
    max_subs = {}
    for i, from_elem in enumerate(atype):
        for j, to_elem in enumerate(atype):
            if i != j:
                max_subs[(from_elem, to_elem)] = min(parent_nat[i], subs_max)

    # ---------- all combinations of substitution counts for each pair from 0 to max_subs
    pairs = list(max_subs.keys())
    max_counts = [max_subs[pair] for pair in pairs]

    # ---------- precompute element index and index pairs for speed
    idx = {elem: i for i, elem in enumerate(atype)}
    index_pairs = [(idx[from_elem], idx[to_elem]) for (from_elem, to_elem) in pairs]

    for counts in product(*[range(max_count + 1) for max_count in max_counts]):
        total_subs = sum(counts)
        if 0 < total_subs <= subs_max:
            # ------ new_nat after substitution
            new_nat = list(parent_nat)
            remain_nat = list(parent_nat)    # to track remaining atoms
            for (from_idx, to_idx), count in zip(index_pairs, counts):
                new_nat[from_idx] -= count
                remain_nat[from_idx] -= count
                if remain_nat[from_idx] < 0:
                    break
                new_nat[to_idx] += count

            # ------ skip if any element goes negative
            if any(n < 0 for n in remain_nat):
                continue

            # ------ check if new_nat is within limits
            if all(ll <= n <= ul for n, ll, ul in zip(new_nat, ll_nat, ul_nat)):
                # ------ charge-neutral check using precomputed cn_comb
                if cn_set is not None:
                    if tuple(new_nat) not in cn_set:
                        continue    # not charge-neutral --> skip this pattern

                # ------ charge-neutral check without precomputed cn_comb
                elif charge is not None:
                    if not is_charge_neutral_nat(new_nat, charge):
                        continue    # not charge-neutral --> skip this pattern

                # ------ OK: store this substitution pattern
                subs_list = []
                for pair, count in zip(pairs, counts):
                    subs_list.extend([pair] * count)
                if subs_list:
                    subs_comb.append(subs_list)

    # ---------- return
    return subs_comb
