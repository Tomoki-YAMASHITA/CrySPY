from dataclasses import dataclass
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
from ...util.struc_util import sort_by_atype, check_distance, get_nat


logger = getLogger('cryspy')


@dataclass(frozen=True)
class PermutationContext:
    """Parameters for permutation generation."""

    atype: tuple[str, ...]
    mindist: tuple[tuple[float, ...], ...]
    ntimes: int
    maxcnt_ea: int


def gen_permutation(
        atype,
        mindist,
        struc_data,
        sp,
        n_perm,
        id_start=None,
        symprec=0.01,
        ntimes=1,
        maxcnt_ea=50,
        rng=None
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_perm (int): number of structures to generate by permutation
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    ntimes (int): number of swaps
    maxcnt_ea (int): maximum number of trial in permutation

    # ---------- return
    offspring_data (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'permutation'}
    '''

    # ---------- initialize
    struc_cnt = 0
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- id_offset
    if id_start is None:
        cid = max(struc_data.keys()) + 1
    else:
        if id_start < (max(struc_data.keys()) + 1):
            logger.error('id_start is already included in structure ID of the data')
        else:
            cid = id_start

    # ---------- generate structures by permutation
    while struc_cnt < n_perm:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ check nat for vc
        if len(parent_A.composition) == 1:
            logger.warning(f'Permutation: {pid_A} is composed of only single element')
            continue
        # ------ generate offspring
        offspring = gen_offspring(
            atype,
            mindist,
            parent_A,
            ntimes,
            maxcnt_ea,
            rng,
        )
        # ------ success
        if offspring is not None:
            offspring_data[cid] = offspring
            parents[cid] = (pid_A, )    # tuple
            operation[cid] = 'permutation'
            try:
                spg_sym, spg_num = offspring.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(offspring, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                    f' from {pid_A:>6} by permutation.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return offspring_data, parents, operation


def gen_offspring(atype, mindist, parent_A, ntimes=1, maxcnt_ea=50, rng=None):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    ntimes (int): number of swaps
    maxcnt_ea (int): maximum number of trial in permutation
    rng (np.random.Generator): random number generator

    # ---------- return
    (if success) offspring (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    parent_species = tuple(site.species_string for site in parent_A)
    cnt = 0

    # ---------- ntimes permutation
    while True:
        offspring = parent_A.copy()    # keep original structure
        n = ntimes
        while n > 0:
            # ------ prepare index for each atom type
            indx_each_type = []
            for a in atype:
                indx_each_type.append(
                    [i for i, site in enumerate(offspring)
                        if site.species_string == a])
            # ------ choose two atom type
            non_empty_indices = [i for i, sublist in enumerate(indx_each_type) if sublist]
            type_choice = rng.choice(non_empty_indices, 2, replace=False)
            # ------ choose index
            indx_choice = []
            for tc in type_choice:
                indx_choice.append(rng.choice(indx_each_type[tc]))
            # ------ replace each other
            offspring.replace(indx_choice[0],
                                species=atype[type_choice[1]])
            offspring.replace(indx_choice[1],
                                species=atype[type_choice[0]])
            n -= 1

        # ------ compare to original one
        offspring_species = tuple(site.species_string for site in offspring)
        if offspring_species == parent_species:
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Permutation: could not change structure' +
                        f' in {maxcnt_ea} times')
                logger.warning('Change parent')
                return None    # change parent
            continue

        # ------ check distance
        success, mindist_ij, dist = check_distance(offspring, atype, mindist)
        if success:
            offspring = sort_by_atype(offspring, atype)
            return offspring
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in permutation: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                logger.warning('Permutation: could not satisfy min_dist' +
                        f' in {maxcnt_ea} times')
                logger.warning('Change parent')
                return None    # change parent


class PermutationOffspringGenerator:
    """Generate one offspring by permutation."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: PermutationContext,
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
            operation='permutation',
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
        parent = self.parent_data.structures[parent_ids[0]]

        # ---------- check nat for vc
        if len(parent.composition) == 1:
            logger.warning(
                f'Permutation: {parent_ids[0]} is composed of only single element'
            )
            return None

        # ---------- generate offspring
        return gen_offspring(
            atype=self.context.atype,
            mindist=self.context.mindist,
            parent_A=parent,
            ntimes=self.context.ntimes,
            maxcnt_ea=self.context.maxcnt_ea,
            rng=rng,
        )


def gen_permutation_batch(
        generator,
        n_perm,
        id_start,
        symprec,
        atype,
        base_seed,
        rng,
    ):
    """Generate permutation offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_perm):
        task_rng = make_task_rng(base_seed, cid, rng)
        result = generator.generate(task_rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring_data[cid] = result.structure
        parents[cid] = result.parent_ids
        operation[cid] = 'permutation'
        try:
            spg_sym, spg_num = result.structure.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(result.structure, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} by permutation.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation
