from dataclasses import dataclass
from logging import getLogger
from typing import Optional

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie

from ..ea_offspring import (
    GenerationResult,
    OffspringResult,
    ParentData,
    SelectionContext,
    generate_with_parent_attempts,
    make_task_rng,
)
from ...util.struc_util import sort_by_atype, get_nat
from ...util.struc_util import check_distance


logger = getLogger('cryspy')


@dataclass(frozen=True)
class StrainContext:
    """Parameters for strain generation."""

    atype: tuple[str, ...]
    mindist: tuple[tuple[float, ...], ...]
    sigma_st: float
    maxcnt_ea: int


def gen_strain(
        atype,
        mindist,
        struc_data,
        sp,
        n_strain,
        id_start=None,
        symprec=0.01,
        sigma_st=0.5,
        maxcnt_ea=50,
        rng=None,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_strain (int): number of structures to generate by strain
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    sigma_st (float): standard deviation for strain matrix
    maxcnt_ea (int): maximum number of trial in strain
    rng (Random Generator): instance of numpy random Generator

    # ---------- return
    offspring_data (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, )}
    operation (dict): {id: 'strain'}
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

    # ---------- generate structures by strain
    while struc_cnt < n_strain:
        # ------ select parents
        pid_A, = sp.get_parents(n_parent=1)    # comma for list[0]
        parent_A = struc_data[pid_A]
        # ------ generate offspring
        offspring = gen_offspring(
            atype,
            mindist,
            parent_A,
            sigma_st,
            maxcnt_ea,
            rng,
        )
        # ------ success
        if offspring is not None:
            offspring_data[cid] = offspring
            parents[cid] = (pid_A, )    # tuple
            operation[cid] = 'strain'
            try:
                spg_sym, spg_num = offspring.get_space_group_info(
                    symprec=symprec
                )
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(offspring, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                    f' from {pid_A:>6} by strain.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return offspring_data, parents, operation


def gen_offspring(
        atype,
        mindist,
        parent_A,
        sigma_st=0.5,
        maxcnt_ea=50,
        rng=None,
    ):
    '''

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    parent_A (Structure): pymatgen Structure object
    sigma_st (float): standard deviation for strain matrix
    maxcnt_ea (int): maximum number of trial in strain
    rng (Random Generator): instance of numpy random Generator

    # ---------- return
    (if success) offspring (Structure): pymatgen Structure object
    (if fail) None
    '''
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    offspring = parent_A.copy()    # keep original structure
    lat_mat = offspring.lattice.matrix.T    # lattice vector as matrix
    cnt = 0

    # ---------- generate strained structure
    while True:
        # ------ strain matrix
        strain_matrix = np.empty([3, 3])
        for i in range(3):
            for j in range(3):
                if i <= j:
                    if i == j:
                        strain_matrix[i][j] = 1.0 + rng.normal(loc=0.0, scale=sigma_st)
                    else:
                        strain_matrix[i][j] = rng.normal(loc=0.0, scale=sigma_st)/2.0
                        strain_matrix[j][i] = strain_matrix[i][j]
        # ------ strained lattice
        strained_lattice = np.dot(strain_matrix, lat_mat).T
        # ------ offspring
        offspring = Structure(strained_lattice, offspring.species, offspring.frac_coords)
        # ------ scale lattice
        offspring.scale_lattice(parent_A.volume)
        # ------ check distance
        success, mindist_ij, dist = check_distance(offspring, atype, mindist)
        if success:
            offspring = sort_by_atype(offspring, atype)
            return offspring
        else:
            type0 = atype[mindist_ij[0]]
            type1 = atype[mindist_ij[1]]
            logger.warning(f'mindist in strain: {type0} - {type1}, {dist}. retry.')
            cnt += 1
            if cnt >= maxcnt_ea:
                return None    # change parent


class StrainOffspringGenerator:
    """Generate one offspring by strain."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: StrainContext,
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
            operation='strain',
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

        # ---------- generate offspring
        return gen_offspring(
            atype=self.context.atype,
            mindist=self.context.mindist,
            parent_A=parent,
            sigma_st=self.context.sigma_st,
            maxcnt_ea=self.context.maxcnt_ea,
            rng=rng,
        )


def gen_strain_batch(
        generator,
        n_strain,
        id_start,
        symprec,
        atype,
        base_seed,
        rng,
    ):
    """Generate strain offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_strain):
        task_rng = make_task_rng(base_seed, cid, rng)
        result = generator.generate(task_rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring_data[cid] = result.structure
        parents[cid] = result.parent_ids
        operation[cid] = 'strain'
        try:
            spg_sym, spg_num = result.structure.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(result.structure, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} by strain.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation