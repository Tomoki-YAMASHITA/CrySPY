from logging import getLogger
from logging.handlers import QueueListener
import multiprocessing as mp
from queue import Empty

import numpy as np
from pymatgen.core import Structure

from ...EA.ea_offspring import (
    EAOffspringGenerator,
    OffspringTask,
    ParentData,
    TournamentSelectionContext,
    build_roulette_selection_context,
)
from ...EA.gen_struc_EA.crossover import (
    CrossoverContext,
    CrossoverOffspringGenerator,
)
from ...EA.gen_struc_EA.permutation import (
    PermutationContext,
    PermutationOffspringGenerator,
)
from ...EA.gen_struc_EA.strain import (
    StrainContext,
    StrainOffspringGenerator,
)
from ...EA.natural_selection import natural_selection
from ...util.struc_util import set_mindist
from .controller_rs import generate_structures_rs
from .worker_ea import run_worker_ea
from ..db.convert import raw_to_struc
from ..db.ea import (
    insert_ea_elites,
    insert_ea_info,
    insert_ea_origin,
    select_ea_elite_ids,
)
from ..db.record import (
    Status,
    get_next_record_id,
    initialize_record,
    select_ids_by_generation,
    select_results,
    select_status,
    select_struc,
    update_generation,
    update_init_struc,
)
from ..db.sqlite import connect_db


logger = getLogger('cryspy')


def load_current_generation_ids(
    conn,
    generation: int,
) -> list[int]:
    """Load structure IDs in current EA generation."""

    # ---------- current generation IDs
    current_ids = select_ids_by_generation(conn, generation)

    # ---------- check IDs
    if not current_ids:
        logger.error(
            'No structures were found for EA generation '
            f'{generation}'
        )
        raise SystemExit(1)

    # ---------- return
    return current_ids


def load_generation_results(
    conn,
    record_ids: list[int],
) -> tuple[dict[int, Structure], dict[int, float]]:
    """Load optimized structures and fitness values."""

    # ---------- load generation results
    results = select_results(
        conn,
        ids=record_ids,
        sort_by='id',
    )
    struc_data = {}
    fitness = {}

    # ---------- load optimized structures
    for (
        cid,
        status,
        energy,
        _init_spg_num,
        _init_spg_sym,
        _opt_spg_num,
        _opt_spg_sym,
        _nat,
    ) in results:
        if status not in (Status.DONE, Status.NOT_CONVERGED):
            continue
        if energy is None:
            continue

        raw_struc = select_struc(conn, cid, initial=False)
        if raw_struc is None:
            continue
        if raw_struc['lattice'] is None:
            continue
        if raw_struc['frac_coords'] is None:
            continue

        struc_data[cid] = raw_to_struc(raw_struc)
        fitness[cid] = energy

    # ---------- check generation results
    if not struc_data:
        logger.error('No generation results were found for EA generation')
        raise SystemExit(1)
    if set(struc_data) != set(record_ids):
        missing_ids = sorted(set(record_ids) - set(struc_data))
        logger.error(
            'Failed to load generation results: '
            f'missing IDs = {missing_ids}'
        )
        raise SystemExit(1)

    # ---------- return
    return struc_data, fitness


def prepare_parent_pool(
    conn,
    rin,
    generation: int,
    rng: np.random.Generator,
):
    """Prepare parent pool for fixed-composition EA."""

    # ---------- current generation
    current_ids = load_current_generation_ids(conn, generation)
    struc_data, fitness = load_generation_results(conn, current_ids)

    # ---------- elite
    elite_ids = select_ea_elite_ids(conn, generation)
    if elite_ids:
        elite_struc, elite_fitness = load_generation_results(
            conn,
            elite_ids,
        )
    else:
        elite_struc = None
        elite_fitness = None

    # ---------- natural selection
    ranking, fit_with_elite, struc_with_elite = natural_selection(
        fitness=fitness,
        struc_data=struc_data,
        elite_struc=elite_struc,
        elite_fitness=elite_fitness,
        n_fittest=rin.n_fittest,
        fit_reverse=rin.fit_reverse,
        emax_ea=rin.emax_ea,
        emin_ea=rin.emin_ea,
        rng=rng,
    )

    # ---------- check parent pool
    min_required = 2 if rin.n_crsov > 0 else 1
    if len(ranking) < min_required:
        logger.error(
            'Not enough parent candidates after natural selection: '
            f'required {min_required}, got {len(ranking)}'
        )
        raise SystemExit(1)

    # ---------- return
    return ranking, fit_with_elite, struc_with_elite


def prepare_parent_data(
    struc_data: dict[int, Structure],
) -> ParentData:
    """Prepare parent structure data."""

    # ---------- parent data
    return ParentData(
        structures=struc_data,
    )


def prepare_selection_context(
    rin,
    ranking: list[int],
    fitness: dict[int, float],
):
    """Prepare parent selection context."""

    # ---------- tournament
    if rin.slct_func == 'TNM':
        return TournamentSelectionContext(
            parent_ids=tuple(ranking),
            tournament_size=rin.t_size,
        )

    return build_roulette_selection_context(
        parent_ids=tuple(ranking),
        fitness=fitness,
        a_rlt=rin.a_rlt,
        b_rlt=rin.b_rlt,
        fit_reverse=rin.fit_reverse,
    )


def prepare_crossover_context(
    rin,
    mindist,
) -> CrossoverContext:
    """Prepare crossover context."""

    # ---------- crossover context
    return CrossoverContext(
        atype=rin.atype,
        nat=rin.nat,
        mindist=mindist,
        crs_lat=rin.crs_lat,
        nat_diff_tole=rin.nat_diff_tole,
        maxcnt_ea=rin.maxcnt_ea,
        vc=False,
        ll_nat=None,
        ul_nat=None,
        cn_comb=None,
        feasible_N=None,
        charge=None,
        cn_data=None,
        min_comp=None,
        max_comp=None,
    )


def prepare_permutation_context(
    rin,
    mindist,
) -> PermutationContext:
    """Prepare permutation context."""

    # ---------- permutation context
    return PermutationContext(
        atype=rin.atype,
        mindist=mindist,
        ntimes=rin.ntimes,
        maxcnt_ea=rin.maxcnt_ea,
    )


def prepare_strain_context(
    rin,
    mindist,
) -> StrainContext:
    """Prepare strain context."""

    # ---------- strain context
    return StrainContext(
        atype=rin.atype,
        mindist=mindist,
        sigma_st=rin.sigma_st,
        maxcnt_ea=rin.maxcnt_ea,
    )


def create_offspring_tasks(
    conn,
    rin,
    generation: int,
) -> list[OffspringTask]:
    """Create offspring generation tasks."""

    # ---------- initialize
    tasks = []
    cid = get_next_record_id(conn)

    # ---------- create tasks
    for operation, count in (
        ('crossover', rin.n_crsov),
        ('permutation', rin.n_perm),
        ('strain', rin.n_strain),
    ):
        for _ in range(count):
            initialize_record(
                conn,
                cid,
                Status.GENERATING,
                generation=generation,
            )
            tasks.append(
                OffspringTask(
                    cid=cid,
                    generation=generation,
                    operation=operation,
                )
            )
            cid += 1

    # ---------- return
    return tasks


def prepare_random_records_ea(
    conn,
    random_ids: list[int],
    generation: int,
) -> list[int]:
    """Prepare records for EA random structure generation."""

    # ---------- initialize
    generate_ids = []

    # ---------- prepare records
    for cid in random_ids:
        status = select_status(conn, cid)

        # ------ existing record
        if status is not None:
            if status == Status.GENERATING:
                generate_ids.append(cid)
            continue

        # ------ initialize missing record
        initialize_record(
            conn,
            cid,
            Status.GENERATING,
            generation=generation,
        )
        generate_ids.append(cid)

    # ---------- return
    return generate_ids


def prepare_offspring_generator(
    rin,
    parent_data: ParentData,
    selection_context,
    crossover_context: CrossoverContext,
    permutation_context: PermutationContext,
    strain_context: StrainContext,
) -> EAOffspringGenerator:
    """Prepare EA offspring generator."""

    # ---------- operation generators
    return EAOffspringGenerator(
        operation_generators={
            'crossover': CrossoverOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context,
                context=crossover_context,
                max_parent_attempts=rin.max_parent_attempts,
            ).generate,
            'permutation': PermutationOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context,
                context=permutation_context,
                max_parent_attempts=rin.max_parent_attempts,
            ).generate,
            'strain': StrainOffspringGenerator(
                parent_data=parent_data,
                selection_context=selection_context,
                context=strain_context,
                max_parent_attempts=rin.max_parent_attempts,
            ).generate,
        },
    )


def run_controller_ea(
    tasks: list[OffspringTask],
    task_queue,
    result_queue,
    stop_event,
    processes: list[mp.Process],
):
    """Control EA offspring generation."""

    # ---------- initialize
    num_active = 0
    worker_failed = False
    results = []

    # ---------- queue tasks
    for task in tasks:
        task_queue.put(task)
        num_active += 1

    # ---------- process results
    while num_active > 0 or worker_failed:
        result = None
        try:
            # ------ get result
            result = result_queue.get(timeout=1.0)

        except Empty:
            pass

        # ------ check workers
        if not worker_failed:
            stopped_processes = [
                process for process in processes
                if process.exitcode is not None
            ]
            if stopped_processes:
                for process in stopped_processes:
                    logger.error(
                        f'{process.name} stopped: '
                        f'exitcode = {process.exitcode}'
                    )

                # ------ drain workers
                worker_failed = True
                stop_event.set()
                for _ in processes:
                    task_queue.put(None)

        # ------ no result
        if result is None:
            if (
                worker_failed
                and all(
                    process.exitcode is not None
                    for process in processes
                )
            ):
                break
            continue

        task, offspring_result, error_message = result
        num_active -= 1

        # ------ check generation error
        if error_message is not None:
            raise RuntimeError(
                f'Offspring generation failed: '
                f'ID {task.cid}: {error_message}'
            )

        # ------ store result
        results.append((task, offspring_result))

    # ---------- worker failure
    if worker_failed:
        raise RuntimeError(
            'Worker stopped before returning all offspring'
        ) from None

    # ---------- return
    return results


def save_offspring_results(
    conn,
    rin,
    results,
) -> None:
    """Save EA offspring generation results."""

    # ---------- save results
    for task, offspring_result in results:
        parent_ids = offspring_result.parent_ids
        parent1_id = parent_ids[0] if len(parent_ids) > 0 else None
        parent2_id = parent_ids[1] if len(parent_ids) > 1 else None

        update_init_struc(
            conn,
            task.cid,
            offspring_result.structure,
            rin.atype,
            rin.symprec,
            Status.WAITING,
        )
        insert_ea_origin(
            conn,
            task.cid,
            task.generation,
            task.operation,
            parent1_id,
            parent2_id,
            offspring_result.parent_attempt_count,
        )


def save_ea_generation_info(
    conn,
    rin,
    generation: int,
) -> None:
    """Save EA generation information."""

    # ---------- save generation information
    insert_ea_info(
        conn,
        generation,
        rin.n_pop,
        rin.n_crsov,
        rin.n_perm,
        rin.n_strain,
        rin.n_rand,
        rin.n_elite,
        rin.crs_lat,
        rin.slct_func,
    )


def save_ea_elite_ids(
    conn,
    generation: int,
    elite_ids: list[int],
) -> None:
    """Save EA elite IDs."""

    # ---------- no elite
    if not elite_ids:
        return

    # ---------- save elite IDs
    insert_ea_elites(
        conn,
        generation,
        elite_ids,
    )


def select_elite_ids(
    rin,
    struc_data: dict[int, Structure],
    fitness: dict[int, float],
    rng: np.random.Generator,
) -> list[int]:
    """Select elite IDs for next generation."""

    # ---------- no elite
    if rin.n_elite <= 0:
        return []

    # ---------- select elite
    ranking, _, _ = natural_selection(
        fitness=fitness,
        struc_data=struc_data,
        elite_struc=None,
        elite_fitness=None,
        n_fittest=rin.n_elite,
        fit_reverse=rin.fit_reverse,
        emax_ea=rin.emax_ea,
        emin_ea=rin.emin_ea,
        rng=rng,
    )

    # ---------- log elite
    for cid in ranking:
        logger.info(f'Structure ID {cid:>6} keeps as the elite')

    # ---------- return
    return ranking


def generate_offspring_ea(
    rin,
    tasks: list[OffspringTask],
    offspring_generator: EAOffspringGenerator,
):
    """Generate EA offspring with workers."""

    # ---------- nothing to generate
    if not tasks:
        return []

    # ---------- number of workers
    num_workers = min(rin.njob, len(tasks))
    logger.info(f'Number of EA generation workers: {num_workers}')

    # ---------- queues and event
    task_queue = mp.Queue()
    result_queue = mp.Queue()
    log_queue = mp.Queue()
    stop_event = mp.Event()

    # ---------- logging
    listener = QueueListener(
        log_queue,
        *logger.handlers,
        respect_handler_level=True,
    )
    listener_started = False

    # ---------- initialize workers
    processes = []

    try:
        # ---------- start workers
        for worker_id in range(num_workers):
            process = mp.Process(
                target=run_worker_ea,
                args=(
                    offspring_generator,
                    rin.seed,
                    task_queue,
                    result_queue,
                    stop_event,
                    log_queue,
                    logger.level,
                ),
                name=f'worker-ea-{worker_id + 1}',
            )
            process.start()
            processes.append(process)
            logger.info(
                f'Launched {process.name}: PID = {process.pid}'
            )

        # ---------- start listener
        listener.start()
        listener_started = True

        # ---------- run controller
        results = run_controller_ea(
            tasks,
            task_queue,
            result_queue,
            stop_event,
            processes,
        )

        # ---------- stop workers
        for _ in processes:
            task_queue.put(None)

        # ---------- wait workers
        for process in processes:
            process.join()

    except BaseException as error:
        # ---------- terminate workers
        for process in processes:
            if process.is_alive():
                process.terminate()
        for process in processes:
            process.join()

        # ---------- execution error
        if isinstance(error, Exception):
            logger.exception(
                'EA offspring generation failed'
            )
            raise SystemExit(1)

        raise

    finally:
        # ---------- stop logging
        if listener_started:
            listener.stop()

        # ---------- close queues
        task_queue.close()
        result_queue.close()
        log_queue.close()
        task_queue.join_thread()
        result_queue.join_thread()
        log_queue.join_thread()

    # ---------- check exit codes
    failed_processes = [
        process for process in processes
        if process.exitcode != 0
    ]
    if failed_processes:
        for process in failed_processes:
            logger.error(
                f'{process.name} failed: '
                f'exitcode = {process.exitcode}'
            )
        raise SystemExit(1)

    # ---------- return
    return results


def generate_and_save_offspring_ea(
    conn,
    rin,
    tasks: list[OffspringTask],
    offspring_generator: EAOffspringGenerator,
    current_struc_data: dict[int, Structure],
    current_fitness: dict[int, float],
    rng: np.random.Generator,
) -> list[int]:
    """Generate and save EA offspring."""

    # ---------- generate offspring
    results = generate_offspring_ea(
        rin,
        tasks,
        offspring_generator,
    )

    # ---------- select elite
    elite_ids = select_elite_ids(
        rin,
        current_struc_data,
        current_fitness,
        rng,
    )

    try:
        # ---------- save results
        save_offspring_results(
            conn,
            rin,
            results,
        )

        conn.commit()

    except Exception:
        # ---------- rollback
        conn.rollback()
        raise

    # ---------- return
    return elite_ids


def initialize_ea(
    rin,
    cids: list[int],
) -> None:
    """Initialize fixed-composition EA."""

    # ---------- generation
    generation = 1

    # ---------- initialize database connection
    conn = None
    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- prepare random records
        generate_ids = prepare_random_records_ea(
            conn,
            cids,
            generation,
        )
        conn.commit()

        # ---------- generate random structures
        generate_structures_rs(
            rin,
            generate_ids,
        )


        # ---------- update generation
        update_generation(
            conn,
            cids,
            generation,
        )

        # ---------- check generation IDs
        generation_ids = set(
            select_ids_by_generation(
                conn,
                generation,
            )
        )
        expected_ids = set(cids)
        if generation_ids != expected_ids:
            missing_ids = sorted(expected_ids - generation_ids)
            unexpected_ids = sorted(generation_ids - expected_ids)
            logger.error(
                f'Initial EA population does not match n_pop: '
                f'missing IDs = {missing_ids}, '
                f'unexpected IDs = {unexpected_ids}'
            )
            raise SystemExit(1)

        # ---------- save random origins
        for cid in cids:
            insert_ea_origin(
                conn,
                cid,
                generation,
                'random',
            )

        # ---------- save generation information
        insert_ea_info(
            conn,
            generation,
            rin.n_pop,
            0,
            0,
            0,
            rin.n_pop,
            0,
            rin.crs_lat,
            rin.slct_func,
        )

        conn.commit()
    except Exception:
        if conn is not None:
            conn.rollback()
        raise
    finally:
        if conn is not None:
            conn.close()


def generate_next_generation_ea(
    rin,
    current_generation: int,
) -> None:
    """Generate and save the next EA generation."""

    # ---------- generation
    next_generation = current_generation + 1
    logger.info('# ---------- Evolutionary algorithm')
    logger.info(f'Generation {next_generation}')

    # ---------- initialize database connection
    conn = None

    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- controller RNG
        if rin.seed is None:
            rng = np.random.default_rng()
        else:
            seed_sequence = np.random.SeedSequence(
                [
                    rin.seed,
                    current_generation,
                ]
            )
            rng = np.random.default_rng(seed_sequence)

        # ---------- parent pool
        ranking, fitness, struc_data = prepare_parent_pool(
            conn,
            rin,
            current_generation,
            rng,
        )

        # ---------- mindist
        logger.info('# ------ mindist')
        mindist = set_mindist(
            rin.atype,
            rin.mindist,
            rin.mindist_factor,
            rin.struc_mode,
        )

        # ---------- generation contexts
        parent_data = prepare_parent_data(struc_data)
        selection_context = prepare_selection_context(
            rin,
            ranking,
            fitness,
        )
        crossover_context = prepare_crossover_context(
            rin,
            mindist,
        )
        permutation_context = prepare_permutation_context(
            rin,
            mindist,
        )
        strain_context = prepare_strain_context(
            rin,
            mindist,
        )

        # ---------- offspring tasks
        tasks = create_offspring_tasks(
            conn,
            rin,
            next_generation,
        )

        # ---------- random IDs
        random_id_start = get_next_record_id(conn)
        random_ids = list(
            range(
                random_id_start,
                random_id_start + rin.n_rand,
            )
        )

        # ---------- population
        if len(tasks) + len(random_ids) != rin.n_pop:
            logger.error(
                'Number of offspring tasks and random IDs '
                f'does not match n_pop: '
                f'{len(tasks)} + {len(random_ids)} '
                f'!= {rin.n_pop}'
            )
            raise SystemExit(1)

        # ---------- offspring generator
        offspring_generator = prepare_offspring_generator(
            rin,
            parent_data,
            selection_context,
            crossover_context,
            permutation_context,
            strain_context,
        )

        # ---------- generate and save offspring
        elite_ids = generate_and_save_offspring_ea(
            conn,
            rin,
            tasks,
            offspring_generator,
            struc_data,
            fitness,
            rng,
        )

        # ---------- prepare random records
        generate_ids = prepare_random_records_ea(
            conn,
            random_ids,
            next_generation,
        )
        conn.commit()

        # ---------- generate random structures
        generate_structures_rs(
            rin,
            generate_ids,
        )

        try:
            # ---------- update random records
            update_generation(
                conn,
                random_ids,
                next_generation,
            )

            # ---------- save random origins
            for cid in random_ids:
                insert_ea_origin(
                    conn,
                    cid,
                    next_generation,
                    'random',
                )

            # ---------- save elite IDs
            save_ea_elite_ids(
                conn,
                next_generation,
                elite_ids,
            )

            # ---------- save generation information
            save_ea_generation_info(
                conn,
                rin,
                next_generation,
            )

            conn.commit()

        except Exception:
            # ---------- rollback
            conn.rollback()
            raise

    finally:
        # ---------- close database
        if conn is not None:
            conn.close()