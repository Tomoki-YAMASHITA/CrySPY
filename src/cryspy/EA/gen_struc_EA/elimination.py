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
from ...util.struc_util import get_nat

# ---------- import later
#from ...util.charge_neutral import is_charge_neutral_nat


logger = getLogger('cryspy')


@dataclass(frozen=True)
class EliminationContext:
    """Parameters for elimination generation."""

    atype: tuple[str, ...]
    elim_dnat_map: dict[int, list[tuple[int, ...]]]


def gen_child(parent_A, elim_element_list, rng=None):
    '''
    parent_A (Structure): pymatgen Structure object
    elim_element_list (list): list of atom types to eliminate, e.g. ['Li', 'Li', 'O']
    rng (np.random.Generator): random number generator

    # ---------- return
    child (Structure): pymatgen Structure object
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- keep original structure
    child = parent_A.copy()    # keep original structure

    # ---------- elimination
    # ------ get indices of atoms to eliminate
    elim_indices = []
    used = set()
    for elem in elim_element_list:
        # まだ使われていないインデックスだけを候補に
        candidates = [i for i, site in enumerate(child) if site.species_string == elem and i not in used]
        idx = rng.choice(candidates)
        elim_indices.append(idx)
        used.add(idx)

    # ------ remove atom
    child.remove_sites(elim_indices)

    # ---------- return
    #            no need to check distance in elimination
    return child


class EliminationOffspringGenerator:
    """Generate one offspring by elimination."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: EliminationContext,
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
            operation='elimination',
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

        # ---------- delta nat combinations
        dnat_comb = self.context.elim_dnat_map[parent_id]
        dnat = dnat_comb[rng.integers(len(dnat_comb))]
        elim_element_list = [
            a for a, n in zip(self.context.atype, dnat)
            for _ in range(n)
        ]

        # ---------- generate offspring
        return gen_child(
            parent_A=parent,
            elim_element_list=elim_element_list,
            rng=rng,
        )


def gen_elimination_batch(
        generator,
        n_elim,
        id_start,
        symprec,
        atype,
        base_seed,
        rng,
    ):
    """Generate elimination offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_elim):
        task_rng = make_task_rng(base_seed, cid, rng)
        result = generator.generate(task_rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring_data[cid] = result.structure
        parents[cid] = result.parent_ids
        operation[cid] = 'elimination'
        try:
            spg_sym, spg_num = result.structure.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(result.structure, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} by elimination.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation


def get_elim_dnat_comb(ll_nat, elim_max, parent_nat, cn_comb=None, charge=None):
    """
    Generate candidate delta nat (number of atoms) for elimination.

    This function returns all feasible atom-count decrements to be eliminated
    from `parent_nat` under the lower-bound and maximum-elimination constraints.

    Parameters
    ----------
    ll_nat : tuple[int]
        Per-element lower bounds of atom counts.
    elim_max : int
        Maximum total number of atoms to eliminate in one operation.
    parent_nat : tuple[int]
        nat of the parent structure.
    cn_comb (np.ndarray): charge neutral combinations
    charge (tuple): charge of each atom type

    Returns
    -------
    dnat_comb : list[tuple[int, ...]]
        Each candidate satisfies:
        - 0 < sum(dnat) <= elim_max
        - parent_nat - dnat >= ll_nat (element-wise)
        - sum(parent_nat - dnat) > 0
        - parent_nat - dnat is charge-neutral if cn_comb or charge is given

    Case 1: `cn_comb is None`
    All dnat candidates are locally enumerated under elim_max and ll_nat.

    Case 2: `cn_comb` is given
    All dnat candidates are locally enumerated, and only candidates where
    parent_nat - dnat is in cn_comb are kept.

    Case 3: `cn_comb is None` and `charge` is given
    All dnat candidates are locally enumerated, and only candidates where
    parent_nat - dnat is charge-neutral are kept.
    """
    # ---------- charge-neutral checker
    if charge is not None and cn_comb is None:
        from ...util.charge_neutral import is_charge_neutral_nat

    # ---------- initialize
    dnat_comb = []
    cn_set = {tuple(row) for row in cn_comb} if cn_comb is not None else None

    # ---------- enumerate elimination candidates
    max_del_per_element = [
        min(current - ll, elim_max)
        for current, ll in zip(parent_nat, ll_nat)
    ]
    for dnat in product(*[range(max_del + 1) for max_del in max_del_per_element]):
        if 0 < sum(dnat) <= elim_max:
            new_nat = tuple(a - b for a, b in zip(parent_nat, dnat))
            if sum(new_nat) <= 0:
                continue

            # ------ charge-neutral check
            if cn_set is not None:
                if new_nat not in cn_set:
                    continue
            elif charge is not None:
                if not is_charge_neutral_nat(new_nat, charge):
                    continue

            dnat_comb.append(dnat)

    # ---------- return
    return dnat_comb
