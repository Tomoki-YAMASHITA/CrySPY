from collections import Counter
from dataclasses import dataclass
from logging import getLogger
from typing import Optional

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.periodic_table import DummySpecie

from ..ea_offspring import (
    GenerationResult,
    OffspringResult,
    ParentData,
    SelectionContext,
    generate_with_parent_attempts,
)
from ...util.struc_util import origin_shift, sort_by_atype, check_distance
from ...util.struc_util import get_nat

# ---------- import later
#from ...util.charge_neutral import is_charge_neutral_nat


logger = getLogger('cryspy')


@dataclass(frozen=True)
class CrossoverContext:
    """Parameters for crossover generation."""

    atype: tuple[str, ...]
    nat: tuple[int, ...] | None
    mindist: tuple[tuple[float, ...], ...]
    crs_lat: str
    nat_diff_tole: int
    maxcnt_ea: int
    vc: bool
    ll_nat: tuple[int, ...] | None
    ul_nat: tuple[int, ...] | None
    cn_comb: np.ndarray | None
    feasible_N: list | None
    charge: tuple[float, ...] | None
    cn_data: dict | None
    min_comp: tuple[float, ...] | None
    max_comp: tuple[float, ...] | None


def gen_crossover(
        atype,
        nat,
        mindist,
        struc_data,
        sp,
        n_crsov,
        id_start=None,
        symprec=0.01,
        crs_lat='random',
        nat_diff_tole=4,
        maxcnt_ea=50,
        vc=False,
        ll_nat=None,
        ul_nat=None,
        cn_comb=None,
        feasible_N=None,
        charge=None,
        cn_data=None,
        min_comp=None,
        max_comp=None,
        rng=None,
    ):
    '''
    # ---------- args

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    nat (tuple): e.g. (4, 4), None if algo == 'EA-vc'
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    struc_data (dict): {id: structure data}
    sp (instance): instance of SelectParents class
    n_crsov (int): number of structures to generate by crossover
    id_start (int): start ID for new structures
    symprec (float): tolerance for symmetry
    crs_lat (str): 'equal' or 'random'
    nat_diff_tole (int): tolerance for nat_diff
    maxcnt_ea (int): maximum number of trial in crossover
    vc (bool): set True if algo == 'EA-vc'
    ll_nat (tuple): lower limit of nat for EA-vc, e.g. (0, 0)
    ul_nat (tuple): upper limit of nat for EA-vc, e.g. (8, 8)
    cn_comb (numpy.ndarray): charge neutral combinations
    feasible_N (list): list of feasible total atom counts for EA-vc with composition constraints
    charge (tuple): charge of each atom type
    cn_data (dict): charge-neutral data for enumerate/sample mode
    min_comp (tuple): lower composition bounds
    max_comp (tuple): upper composition bounds
    rng (numpy.random.Generator): random number generator

    # ---------- return
    offspring_data (dict): {id: structure data}
    parents (dict): {id: (id of parent_A, id of parent_B)}
    operation (dict): {id: 'crossover'}
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
            raise SystemExit(1)
        else:
            cid = id_start

    # ---------- generate structures by crossover
    while struc_cnt < n_crsov:
        # ------ select parents
        pid_A, pid_B = sp.get_parents(n_parent=2)    # get IDs
        parent_A = struc_data[pid_A]
        parent_B = struc_data[pid_B]
        # ------ generate offspring
        offspring = gen_offspring(
            atype,
            nat,
            mindist,
            parent_A,
            parent_B,
            crs_lat,
            nat_diff_tole,
            maxcnt_ea,
            vc,
            ll_nat,
            ul_nat,
            cn_comb,
            feasible_N,
            charge,
            cn_data,
            min_comp,
            max_comp,
            rng,
        )
        # ------ success
        if offspring is not None:
            # -- Niggli reduction
            try:
                offspring = offspring.get_reduced_structure(reduction_algo="niggli")
            except Exception as e:
                logger.warning(f'Niggli reduction failed: {e}')
            offspring_data[cid] = offspring
            parents[cid] = (pid_A, pid_B)
            operation[cid] = 'crossover'
            try:
                spg_sym, spg_num = offspring.get_space_group_info(symprec=symprec)
            except TypeError:
                spg_num = 0
                spg_sym = None
            tmp_nat = get_nat(offspring, atype)
            logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                    f' from {pid_A:>6} and {pid_B:>6} by crossover.'
                    f' Space group: {spg_num:>3} {spg_sym}')
            cid += 1
            struc_cnt += 1

    # ---------- return
    return offspring_data, parents, operation


def gen_offspring(
        atype,
        nat,
        mindist,
        parent_A,
        parent_B,
        crs_lat='random',
        nat_diff_tole=4,
        maxcnt_ea=50,
        vc=False,
        ll_nat=None,
        ul_nat=None,
        cn_comb=None,
        feasible_N=None,
        charge=None,
        cn_data=None,
        min_comp=None,
        max_comp=None,
        rng=None,
    ):
    '''
    # ---------- args

        tuple may be replaced by list

    atype (tuple): e.g. ('Na', 'Cl')
    nat (tuple): e.g. (4, 4), None if algo == 'EA-vc'
    parent_A (Structure): pymatgen Structure object
    parent_B (Structure): pymatgen Structure object
    mindist (tuple): minimum interatomic distance, e.g. ((1.5, 1.5), (1.5, 1.5))
    crs_lat (str): 'equal' or 'random'
    nat_diff_tole (int): tolerance for nat_diff
    maxcnt_ea (int): maximum number of trial in crossover
    vc (bool): set True if algo == 'EA-vc'
    ll_nat (tuple): lower limit of nat for EA-vc, e.g. (1, 1)
    ul_nat (tuple): upper limit of nat for EA-vc, e.g. (8, 8)
    cn_comb (numpy.ndarray): charge neutral combinations
    feasible_N (list): list of feasible total atom counts for EA-vc with composition constraints
    charge (tuple): charge of each atom type
    cn_data (dict): charge-neutral data for enumerate/sample mode
    min_comp (tuple): lower composition bounds
    max_comp (tuple): upper composition bounds
    rng (numpy.random.Generator): random number generator

    # ---------- return
    (if success) offspring (Structure): pymatgen Structure object
    (if fail) None
    '''

    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- initialize
    parent_A = origin_shift(parent_A, rng)    # origin_shift returns a new Structure object
    parent_B = origin_shift(parent_B, rng)
    count = 0

    # ---------- lattice crossover
    if crs_lat == 'equal':
        w_lat = np.array([1.0, 1.0])
    elif crs_lat == 'random':
        w_lat = rng.choice([0.0, 1.0], size=2, replace=False)
    else:
        logger.error('crs_lat must be equal or random')
    lattice = _lattice_crossover(parent_A, parent_B, w_lat)

    # ---------- generate offspring
    while True:
        count += 1
        # ------ coordinate crossover
        axis, slice_point, species, coords = _one_point_crossover(parent_A, parent_B, rng)
        # ------ offspring structure
        offspring = Structure(lattice, species, coords)
        # ------ check nat_diff
        # -- charge neutrality in sample mode
        if vc and cn_data is not None and cn_data['mode'] == 'sample':
            # -- choose close target_nat by local search around offspring_nat
            offspring_nat = get_nat(offspring, atype)
            target_nat = _sample_close_cn_nat(
                offspring_nat,
                ll_nat,
                ul_nat,
                charge,
                nat_diff_tole=nat_diff_tole,
                min_comp=min_comp,
                max_comp=max_comp,
                rng=rng,
            )

            # -- no close charge-neutral target_nat found
            use_target_nat = True
            if target_nat is None:
                if count > maxcnt_ea:
                    return None
                continue    # slice again

            # -- check difference from target_nat
            nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)
            nat_diff_checked = False

        # -- charge neutrality in enumerate mode
        elif vc and cn_comb is not None:
            # -- choose closest target_nat from precomputed charge-neutral combinations
            target_nat = _get_close_cn_comb(offspring, atype, cn_comb, rng)
            use_target_nat = True
            nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)
            nat_diff_checked = False

        # -- composition constraints only
        elif vc and feasible_N is not None:
            # -- choose close target_nat from feasible total atom counts
            target_nat, nat_diff = _get_close_feasible_nat(offspring, atype, feasible_N, nat_diff_tole, rng)
            use_target_nat = True
            nat_diff_checked = True

            # -- no close composition-valid target_nat found
            if target_nat is None:
                if count > maxcnt_ea:
                    return None
                continue    # slice again

        else:
            target_nat = nat
            use_target_nat = False
            nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)
            nat_diff_checked = False
        if (not nat_diff_checked) and any([abs(n) > nat_diff_tole for n in nat_diff]):
            logger.debug(f'nat_diff = {nat_diff}')
            if count > maxcnt_ea:    # fail
                return None
            continue    # slice again
        # ------ check mindist
        # either tmp_atype or atype is OK in check_distance()
        success, _, _ = check_distance(offspring, atype, mindist, check_all=False)
        # ------ something smaller than mindist
        if not success:
            # -- remove atoms within mindist
            if any([n > 0 for n in nat_diff]):
                offspring = _remove_within_mindist(offspring, atype, mindist, nat_diff)
                if offspring is None:    # fail --> slice again
                    if count > maxcnt_ea:
                        return None
                    continue
            else:    # nothing to remove, nat_diff = [0, 0]
                if count > maxcnt_ea:
                    return None
                continue    # fail --> slice again
        # ------ recheck nat_diff
        # ------ excess of atoms
        nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)    # recheck
        if any([n > 0 for n in nat_diff]):
            offspring = _remove_border_line(offspring, atype, axis,
                                        slice_point, nat_diff)
        # ------ lack of atoms
        nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)    # recheck
        if any([n < 0 for n in nat_diff]):
            offspring = _add_border_line(offspring, atype, mindist, axis, slice_point,
                                        nat_diff, maxcnt_ea, rng)
        # ------ success --> break while loop
        if offspring is not None:
            break
        # ------ fail --> slice again
        else:
            if count > maxcnt_ea:
                return None
            continue

    # ---------- final check for nat
    nat_diff = _get_nat_diff(atype, target_nat, offspring, vc, ll_nat, ul_nat, use_target_nat)
    if not all([n == 0 for n in nat_diff]):
        return None    # failure

    # ---------- sort by atype
    offspring = sort_by_atype(offspring, atype)

    # ---------- return
    return offspring


class CrossoverOffspringGenerator:
    """Generate one offspring by crossover."""

    def __init__(
        self,
        parent_data: ParentData,
        selection_context: SelectionContext,
        context: CrossoverContext,
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
            operation='crossover',
            selection_context=self.selection_context,
            n_parent=2,
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

        # ---------- parents
        parent_A = self.parent_data.structures[parent_ids[0]]
        parent_B = self.parent_data.structures[parent_ids[1]]

        # ---------- generate offspring
        return gen_offspring(
            atype=self.context.atype,
            nat=self.context.nat,
            mindist=self.context.mindist,
            parent_A=parent_A,
            parent_B=parent_B,
            crs_lat=self.context.crs_lat,
            nat_diff_tole=self.context.nat_diff_tole,
            maxcnt_ea=self.context.maxcnt_ea,
            vc=self.context.vc,
            ll_nat=self.context.ll_nat,
            ul_nat=self.context.ul_nat,
            cn_comb=self.context.cn_comb,
            feasible_N=self.context.feasible_N,
            charge=self.context.charge,
            cn_data=self.context.cn_data,
            min_comp=self.context.min_comp,
            max_comp=self.context.max_comp,
            rng=rng,
        )


def gen_crossover_batch(
        generator,
        n_crsov,
        id_start,
        symprec,
        atype,
        rng,
    ):
    """Generate crossover offspring batch."""

    # ---------- initialize
    offspring_data = {}
    parents = {}
    operation = {}

    # ---------- generate offspring
    for cid in range(id_start, id_start + n_crsov):
        result = generator.generate(rng)
        if not isinstance(result, OffspringResult):
            logger.error(result.reason)
            raise SystemExit(1)
        offspring = result.structure

        # -- Niggli reduction
        try:
            offspring = offspring.get_reduced_structure(
                reduction_algo="niggli"
            )
        except Exception as e:
            logger.warning(f'Niggli reduction failed: {e}')

        offspring_data[cid] = offspring
        parents[cid] = result.parent_ids
        operation[cid] = 'crossover'
        try:
            spg_sym, spg_num = offspring.get_space_group_info(
                symprec=symprec
            )
        except TypeError:
            spg_num = 0
            spg_sym = None
        tmp_nat = get_nat(offspring, atype)
        logger.info(f'Structure ID {cid:>6} {tmp_nat} was generated'
                f' from {result.parent_ids[0]:>6} and {result.parent_ids[1]:>6} by crossover.'
                f' Space group: {spg_num:>3} {spg_sym}')

    # ---------- return
    return offspring_data, parents, operation


def _lattice_crossover(parent_A, parent_B, w_lat):
    # ---------- component --> w_lat
    matrix = ((w_lat[0]*parent_A.lattice.matrix
                + w_lat[1]*parent_B.lattice.matrix)
                / w_lat.sum())
    mat_len = np.sqrt((matrix**2).sum(axis=1))
    # ---------- absolute value of vector
    lat_len = ((np.array(parent_A.lattice.abc)*w_lat[0]
                + np.array(parent_B.lattice.abc)*w_lat[1])
                / w_lat.sum())
    # ---------- correction of vector length
    lat_array = np.empty([3, 3])
    for i in range(3):
        lat_array[i] = matrix[i]*lat_len[i]/mat_len[i]
    # ---------- Lattice in pymatgen
    return Lattice(lat_array)


def _one_point_crossover(parent_A, parent_B, rng=None):
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- slice point and axis
    slice_point = np.clip(rng.normal(loc=0.5, scale=0.1), 0.3, 0.7)
    axis = rng.choice([0, 1, 2])

    # ---------- crossover
    species_A = []
    species_B = []
    coords_A = []
    coords_B = []
    for i in range(parent_A.num_sites):
        if parent_A.frac_coords[i, axis] <= slice_point:
            species_A.append(parent_A[i].species_string)
            coords_A.append(parent_A[i].frac_coords)
        else:
            species_B.append(parent_A[i].species_string)
            coords_B.append(parent_A[i].frac_coords)
    for i in range(parent_B.num_sites):
        if parent_B.frac_coords[i, axis] >= slice_point:
            species_A.append(parent_B[i].species_string)
            coords_A.append(parent_B[i].frac_coords)
        else:
            species_B.append(parent_B[i].species_string)
            coords_B.append(parent_B[i].frac_coords)

    # ---------- adopt a structure with more atoms
    if len(species_A) > len(species_B):
        species = species_A
        coords = coords_A
    elif len(species_A) < len(species_B):
        species = species_B
        coords = coords_B
    else:
        if rng.choice([0, 1]):
            species = species_A
            coords = coords_A
        else:
            species = species_B
            coords = coords_B

    # ---------- return
    return axis, slice_point, species, coords


def _get_nat_diff(atype, target_nat, child, vc, ll_nat, ul_nat, use_target_nat):
    '''
    if not vc or use_target_nat:
        child nat - target nat
        e.g.
            target_nat = [4, 4]
            tmp_nat = [3, 5]
            nat_diff = [-1, 1]

    if vc and not use_target_nat:
        difference from the allowed range [ll_nat, ul_nat]
        e.g.
            ll_nat = [4, 4, 4]
            ul_nat = [8, 8, 8]
            tmp_nat = [2, 6, 12]
            nat_diff = [-2, 0, 4]
    '''
    tmp_nat = get_nat(child, atype)
    if not vc or use_target_nat:
        nat_diff = [i - j for i, j in zip(tmp_nat, target_nat)]
    else:
        nat_diff = []
        for i, n in enumerate(tmp_nat):
            if n > ul_nat[i]:
                nat_diff.append(n - ul_nat[i])
            elif n < ll_nat[i]:
                nat_diff.append(n - ll_nat[i])
            else:
                nat_diff.append(0)
    return nat_diff


def _remove_within_mindist(child, atype, mindist, nat_diff):
    '''
    if success: return child
    if fail:    return None
    '''
    for itype in range(len(atype)):
        while nat_diff[itype] > 0:
            # ---------- check dist
            dist_list = check_distance(child, atype, mindist, check_all=True)
            if not dist_list:    # nothing within mindist
                return child
            # ---------- appearance frequency
            ij_within_dist = [isite[0] for isite in dist_list] + [
                jsite[1] for jsite in dist_list]
            site_counter = Counter(ij_within_dist)
            # ---------- get index for removing
            rm_index = None
            # ---- site[0]: index, site[1]: count
            for site in site_counter.most_common():
                if child[site[0]].species_string == atype[itype]:
                    rm_index = site[0]
                    break    # break for loop
            # ---------- remove atom
            if rm_index is None:    # fail
                return None
            else:
                child.remove_sites([rm_index])
                nat_diff[itype] -= 1

    # ---------- final check
    dist_list = check_distance(child, atype, mindist, check_all=True)
    if dist_list:    # still something within mindist
        logger.warning('remove_within_mindist: some atoms within mindist. retry.')
        return None
    else:    # success
        return child


def _remove_border_line(child, atype, axis, slice_point, nat_diff):
    # ---------- rank atoms from border line
    coords_axis = child.frac_coords[:, axis]

    # ---------- boundary --> 0.0, slice_point, 1.0
    near_sp = (slice_point/2.0 < coords_axis) & \
        (coords_axis < (slice_point + 1.0)/2.0)
    near_one = (slice_point + 1.0)/2.0 <= coords_axis

    # ---------- distance from nearest boundary
    coords_diff = np.where(near_sp,
                            abs(coords_axis - slice_point),
                            coords_axis)
    coords_diff = np.where(near_one, 1.0 - coords_diff, coords_diff)
    atom_border_indx = np.argsort(coords_diff)

    # ---------- remove list
    rm_list = []
    for itype, nrm in enumerate(nat_diff):
        rm_list.append([])
        if nrm > 0:
            for ab_indx in atom_border_indx:
                if child[ab_indx].species_string == atype[itype]:
                    rm_list[itype].append(ab_indx)
                if len(rm_list[itype]) == nrm:
                    break

    # ---------- remove
    for each_type in rm_list:
        if each_type:
            child.remove_sites(each_type)

    # ---------- return
    return child


def _add_border_line(child, atype, mindist, axis, slice_point, nat_diff, maxcnt_ea=50, rng=None):
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- add atoms
    for i in range(len(atype)):
        # ------ counter
        cnt = 0

        # ------ add atoms
        while nat_diff[i] < 0:
            cnt += 1
            coords = rng.random(3)
            mean = _mean_choice(child, axis, slice_point, rng)
            coords[axis] = rng.normal(loc=mean, scale=0.08)
            child.append(species=atype[i], coords=coords)
            success, mindist_ij, dist = check_distance(child, atype, mindist)
            if success:
                cnt = 0    # reset
                nat_diff[i] += 1
            else:
                type0 = atype[mindist_ij[0]]
                type1 = atype[mindist_ij[1]]
                logger.warning(f'mindist in _add_border_line: {type0} - {type1}, {dist}. retry.')
                child.pop()    # cancel
            # -- fail
            if cnt == maxcnt_ea:
                return None

    # ---------- return
    return child


def _mean_choice(child, axis, slice_point, rng=None):
    '''
    Which border contains the most atoms?
    '''
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()
    # ---------- count
    n_zero = np.sum(np.abs(child.frac_coords[:, axis] - 0.0)
                    < 0.1)
    n_slice = np.sum(np.abs(child.frac_coords[:, axis]
                            - slice_point) < 0.1)
    if n_zero < n_slice:
        mean = 0.0
    elif n_zero > n_slice:
        mean = slice_point
    else:
        mean = rng.choice([0.0, slice_point])
    return mean


def _get_close_cn_comb(child, atype, cn_comb, rng=None):
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- distance between child and cn_comb
    tmp_nat = get_nat(child, atype)
    distances = np.sum(np.abs(cn_comb - tmp_nat), axis=1)

    # ---------- find the closest combination
    min_distance = np.min(distances)
    min_indices = np.where(distances == min_distance)[0]

    # ---------- randomly select one of the closest combinations
    closest_idx = rng.choice(min_indices)
    cn_target_nat = cn_comb[closest_idx]

    # ---------- return
    return cn_target_nat


def _get_close_feasible_nat(child, atype, feasible_N, nat_diff_tole, rng=None):
    # ---------- initialize rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- child_nat
    child_nat_arr = np.array(get_nat(child, atype), dtype=int)
    N_child = int(child_nat_arr.sum())

    # ---------- sort feasible_N by distance to child_nat sum
    sorted_feasible = sorted(
        feasible_N,
        key=lambda x: abs(int(x[0]) - N_child)
    )

    # ---------- find the closest feasible nat within nat_diff_tole
    for N, lower, upper in sorted_feasible:
        lower = np.array(lower, dtype=int)
        upper = np.array(upper, dtype=int)
        target_nat_arr = _project_to_N_bounds(child_nat_arr, int(N), lower, upper)
        if target_nat_arr is None:
            continue
        nat_diff = child_nat_arr - target_nat_arr
        if all(abs(int(v)) <= nat_diff_tole for v in nat_diff):
            return tuple(int(v) for v in target_nat_arr), [int(v) for v in nat_diff]

    # ---------- if no feasible nat within nat_diff_tole
    return None, None


def _project_to_N_bounds(nat_arr, N, lower, upper):
    """
    Minimize L1 distance to nat_arr under:
        lower_i <= n_i <= upper_i, sum_i n_i = N, n_i integer
    """

    # ---------- clip each species count into [lower, upper]
    nat_proj_arr = np.clip(nat_arr, lower, upper).astype(int)
    s = int(nat_proj_arr.sum())

    # Strategy: distribute additions/removals
    #           so that the result stays close to the original nat_arr

    # ---------- need to increase sum
    if s < N:
        need = N - s
        # adjust N by adding first to species that are farthest below the original nat_arr
        order = np.argsort(-(nat_arr - nat_proj_arr))  # descending deficit
        for i in order:
            cap = int(upper[i] - nat_proj_arr[i])
            if cap <= 0:
                continue
            add_i = min(cap, need)
            nat_proj_arr[i] += add_i
            need -= add_i
            if need == 0:
                break

    # ---------- need to decrease sum
    elif s > N:
        need = s - N
        # adjust N by removing first from species that are farthest above the original nat_arr
        order = np.argsort(-(nat_proj_arr - nat_arr))  # descending excess
        for i in order:
            cap = int(nat_proj_arr[i] - lower[i])
            if cap <= 0:
                continue
            sub_i = min(cap, need)
            nat_proj_arr[i] -= sub_i
            need -= sub_i
            if need == 0:
                break

    # ---------- final check
    s = int(nat_proj_arr.sum())
    if s != N:
        logger.error(f'_project_to_N_bounds failed: sum={s}, N={N}')
        return None

    # ---------- return
    return nat_proj_arr


def _sample_close_cn_nat(
        child_nat,
        ll_nat,
        ul_nat,
        charge,
        nat_diff_tole=4,
        min_comp=None,
        max_comp=None,
        rng=None,
        tol=1e-12,
    ):
    """
    Sample a charge-neutral nat close to child_nat.
    """
    # ---------- rng
    if rng is None:
        rng = np.random.default_rng()

    # ---------- convert to np.array
    child_nat = np.array(child_nat, dtype=int)
    ll = np.array(ll_nat, dtype=int)
    ul = np.array(ul_nat, dtype=int)
    k = len(child_nat)

    # ---------- make dnat grid
    vals = np.arange(-nat_diff_tole, nat_diff_tole + 1, dtype=int)
    mesh = np.meshgrid(*([vals] * k), indexing='ij')
    dnat_grid = np.stack(mesh, axis=-1).reshape(-1, k)

    # ---------- distance from child_nat
    dist = np.abs(dnat_grid).sum(axis=1)

    # ---------- target nat grid
    target_grid = child_nat + dnat_grid

    # ---------- check bounds
    mask = np.all(target_grid >= ll, axis=1)
    mask &= np.all(target_grid <= ul, axis=1)

    # ---------- check total atoms
    total = target_grid.sum(axis=1)
    mask &= total > 0

    # ---------- check composition constraints
    if min_comp is not None or max_comp is not None:
        comp = np.zeros_like(target_grid, dtype=float)
        comp[mask] = target_grid[mask] / total[mask, None]
        if min_comp is not None:
            min_comp = np.asarray(min_comp, dtype=float)
            mask &= np.all(comp >= min_comp - tol, axis=1)
        if max_comp is not None:
            max_comp = np.asarray(max_comp, dtype=float)
            mask &= np.all(comp <= max_comp + tol, axis=1)

    # ---------- check charge neutrality
    if all(isinstance(c, int) for c in charge):
        # ------ single-valence case
        charges = np.array(charge, dtype=int)
        mask &= (target_grid @ charges == 0)
        candidates = target_grid[mask]
        candidate_dist = dist[mask]
    else:
        # ------ multi-valence case
        from ...util.charge_neutral import is_charge_neutral_nat
        candidates = []
        candidate_dist = []
        # -- check each candidate in the extended valence space
        for target_nat, d in zip(target_grid[mask], dist[mask]):
            if is_charge_neutral_nat(target_nat, charge):
                candidates.append(target_nat)
                candidate_dist.append(d)
        if len(candidates) == 0:
            return None
        candidates = np.array(candidates, dtype=int)
        candidate_dist = np.array(candidate_dist, dtype=int)

    # ---------- no candidate
    if len(candidates) == 0:
        return None

    # ---------- choose one of the closest candidates
    min_dist = candidate_dist.min()
    closest = candidates[candidate_dist == min_dist]
    nat = closest[rng.integers(len(closest))]

    # ---------- return
    return tuple(int(n) for n in nat)
