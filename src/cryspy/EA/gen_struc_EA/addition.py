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

# ---------- import later
#from ...util.charge_neutral import is_charge_neutral_nat


logger = getLogger('cryspy')


@dataclass(frozen=True)
class AdditionContext:
    """Parameters for addition generation."""

    atype: tuple[str, ...]
    mindist: tuple[tuple[float, ...], ...]
    add_dnat_map: dict[int, list[tuple[int, ...]]]
    maxcnt_ea: int


def gen_child(
        atype,
        mindist,
        parent_A,
        add_element_list,
        maxcnt_ea=50,
        rng=None,
    ):
    '''
        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    add_element_list (list): list of atom types to add, e.g. ['Li', 'Li', 'O']
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
    vol10per = False
    vol20per = False
    child = parent_A.copy()    # keep original structure

    # ---------- generate child
    for at in add_element_list:
        cnt = 0
        while True:
            coords = rng.random(3) 
            child.append(species=at, coords=coords)
            # ------ check distance
            success, mindist_ij, dist = check_distance(child, atype, mindist)
            if success:
                break
                # child = sort_by_atype(child, atype)
                # return child
            else:
                type0 = atype[mindist_ij[0]]
                type1 = atype[mindist_ij[1]]
                logger.warning(f'mindist in addition: {type0} - {type1}, {dist}. retry.')
                cnt += 1
                child.pop()    # remove last atom
                if cnt >= maxcnt_ea:
                    # ------ volume change
                    if not vol10per:
                        vol10per = True
                        cnt = 0
                        child.scale_lattice(parent_A.volume * 1.1)
                        logger.warning('Addition: increase volume by 10%')
                        continue
                    elif vol10per and not vol20per:
                        vol20per = True
                        cnt = 0
                        child.scale_lattice(parent_A.volume * 1.2)
                        logger.warning('Addition: increase volume by 20%')
                        continue
                    # ------ fail
                    logger.warning('Addition: could not satisfy min_dist' +
                            f' in {maxcnt_ea} times')
                    logger.warning('Change parent')
                    return None    # change parent

    # ---------- complete
    child = sort_by_atype(child, atype)
    return child


class AdditionOffspringGenerator:
    """Generate one offspring by addition."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: AdditionContext,
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
            operation='addition',
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
        dnat_comb = self.context.add_dnat_map[parent_id]
        dnat = dnat_comb[rng.integers(len(dnat_comb))]
        add_element_list = [
            a for a, n in zip(self.context.atype, dnat)
            for _ in range(n)
        ]

        # ---------- generate offspring
        return gen_child(
            atype=self.context.atype,
            mindist=self.context.mindist,
            parent_A=parent,
            add_element_list=add_element_list,
            maxcnt_ea=self.context.maxcnt_ea,
            rng=rng,
        )


def gen_addition_batch(
        generator,
        n_add,
        id_start,
        symprec,
        atype,
        base_seed,
        rng,
    ):
    """Generate addition offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_add):
        task_rng = make_task_rng(base_seed, cid, rng)
        result = generator.generate(task_rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring_data[cid] = result.structure
        parents[cid] = result.parent_ids
        operation[cid] = 'addition'
        try:
            spg_sym, spg_num = result.structure.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(result.structure, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} by addition.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation


def get_add_dnat_comb(ul_nat, add_max, parent_nat, cn_comb=None, charge=None):
    """
    Generate candidate delta nat (number of atoms) for addition.

    This function returns all feasible atom-count increments to be added to
    `parent_nat` under the upper-bound and maximum-addition constraints.

    Parameters
    ----------
    ul_nat : tuple[int]
        Per-element upper bounds of atom counts.
    add_max : int
        Maximum total number of atoms to add in one operation.
    parent_nat : tuple[int]
        nat of the parent structure.
    cn_comb (np.ndarray): charge neutral combinations
    charge (tuple): charge of each atom type

    Returns
    -------
    dnat_comb : list[tuple[int, ...]]
        Each candidate satisfies:
        - 0 < sum(dnat) <= add_max
        - parent_nat + dnat <= ul_nat (element-wise)
        - parent_nat + dnat is charge-neutral if cn_comb or charge is given

    Examples
    --------
    Case 1: `cn_comb is None`
    parent_nat = (2, 4, 3), ul_nat = (4, 6, 5), add_max = 2

    Possible outputs include:
    [
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
        (2, 0, 0), (0, 2, 0), (0, 0, 2),
        (1, 1, 0), (1, 0, 1), (0, 1, 1),
    ]

    For example:
    - dnat = (1, 0, 1) -> new_nat = (3, 4, 4)
    - dnat = (0, 2, 0) -> new_nat = (2, 6, 3)

    Case 2: `cn_comb` is given
    All dnat candidates are locally enumerated, and only candidates where
    parent_nat + dnat is in cn_comb are kept.

    Case 3: `cn_comb is None` and `charge` is given
    All dnat candidates are locally enumerated, and only candidates where
    parent_nat + dnat is charge-neutral are kept.
    """
    # ---------- charge-neutral checker
    if charge is not None and cn_comb is None:
        from ...util.charge_neutral import is_charge_neutral_nat

    # ---------- initialize
    dnat_comb = []
    cn_set = {tuple(row) for row in cn_comb} if cn_comb is not None else None

    # ---------- enumerate addition candidates
    max_add_per_element = [
        min(ul - current, add_max)
        for ul, current in zip(ul_nat, parent_nat)
    ]
    for dnat in product(*[range(max_add + 1) for max_add in max_add_per_element]):
        if 0 < sum(dnat) <= add_max:
            new_nat = tuple(a + b for a, b in zip(parent_nat, dnat))

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
