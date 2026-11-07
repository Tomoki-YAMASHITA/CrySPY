from collections.abc import Sequence
from datetime import timedelta
from logging import getLogger
from logging.handlers import QueueListener
import multiprocessing as mp
from pathlib import Path
from queue import Empty
from time import perf_counter

from ...IO.read_input import ReadInput
from ...util.charge_neutral import prepare_cn_data
from ...util.struc_util import get_feasible_composition, precompute_feasible_N
from ..db.record_db import (
    Status,
    initialize_record,
    select_status,
    update_init_struc,
)
from ..db.sqlite_db import connect_db
from .worker_rs import run_worker_rs


logger = getLogger('cryspy')


def prepare_generation_record(
    conn,
    cid: int,
) -> None:
    """Prepare a record for random structure generation."""

    # ---------- record status
    status = select_status(conn, cid)

    # ---------- initialize missing record
    if status is None:
        initialize_record(
            conn,
            cid,
            Status.GENERATING,
        )
        conn.commit()
        return

    # ---------- resume generating record
    if status != Status.GENERATING:
        raise RuntimeError(
            f'Cannot generate existing record: '
            f'ID {cid}, status = {status.name}'
        )


def run_controller_rs(
    rin: ReadInput,
    cids: Sequence[int],
    task_queue,
    result_queue,
    stop_event,
    processes: list[mp.Process],
) -> bool:
    """Control database access and random structure generation."""

    # ---------- initialize database connection
    conn = None

    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- initial tasks
        cid_iter = iter(cids)
        num_active = 0
        worker_failed = False
        stop_requested = False
        for _ in processes:
            if Path('STOP_CRYSPY_HT').is_file():
                stop_requested = True
                logger.info('Stop file detected: STOP_CRYSPY_HT')
                break

            try:
                cid = next(cid_iter)
            except StopIteration:
                break
            prepare_generation_record(
                conn,
                cid,
            )
            task_queue.put(cid)
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

            cid, struc, error_message = result
            num_active -= 1

            # ------ check generation error
            if error_message is not None:
                raise RuntimeError(
                    f'Structure generation failed: '
                    f'ID {cid}: {error_message}'
                )

            try:
                # ------ update initial structure
                update_init_struc(
                    conn,
                    cid,
                    struc,
                    rin.atype,
                    rin.symprec,
                    Status.WAITING,
                    vc=rin.algo == 'EA-vc',
                )
                conn.commit()

            except Exception:
                # ------ rollback
                conn.rollback()
                raise

            # ------ stop assigning tasks
            if worker_failed:
                continue

            # ------ check stop file
            if not stop_requested and Path('STOP_CRYSPY_HT').is_file():
                stop_requested = True
                logger.info('Stop file detected: STOP_CRYSPY_HT')
            if stop_requested:
                continue

            # ------ queue next task
            try:
                next_cid = next(cid_iter)
            except StopIteration:
                continue

            prepare_generation_record(
                conn,
                next_cid,
            )
            task_queue.put(next_cid)
            num_active += 1

        # ---------- worker failure
        if worker_failed:
            raise RuntimeError(
                'Worker stopped before returning all structures'
            ) from None

        # ---------- return
        return stop_requested

    finally:
        # ---------- close database
        if conn is not None:
            conn.close()


def generate_structures_rs(
    rin: ReadInput,
    cids: Sequence[int],
    feasible_N: list | None = None,
    cn_data: dict | None = None,
) -> None:
    """Generate and insert random structures."""

    # ---------- nothing to generate
    if not cids:
        return

    # ---------- check stop file
    if Path('STOP_CRYSPY_HT').is_file():
        logger.info('Stop file detected: STOP_CRYSPY_HT')
        raise SystemExit()

    # ---------- log
    logger.info('# ------ Random structure generation')

    # ---------- time
    if rin.algo == 'RS':
        time_start = perf_counter()

    # ---------- charge neutrality
    if rin.algo == 'EA-vc' and rin.charge is not None:
        if cn_data is None:
            cn_data = prepare_cn_data(
                rin.ll_nat,
                rin.ul_nat,
                rin.charge,
                cn_mode=rin.cn_mode,
                max_cn_grid_points=rin.max_cn_grid_points,
                min_comp=rin.min_comp,
                max_comp=rin.max_comp,
            )
    else:
        cn_data = None

    # ---------- composition constraints
    if feasible_N is None and rin.algo == 'EA-vc' and (
        rin.min_comp is not None or rin.max_comp is not None
    ):
        feasible_comp = get_feasible_composition(rin.min_comp, rin.max_comp)
        feasible_N = precompute_feasible_N(rin.ll_nat, rin.ul_nat, feasible_comp)

    # ---------- number of workers
    njob_struc = rin.njob if rin.njob_struc is None else rin.njob_struc
    num_workers = min(njob_struc, len(cids))
    logger.info(f'Number of generation workers: {num_workers}')

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
                target=run_worker_rs,
                args=(
                    rin,
                    worker_id == 0,
                    task_queue,
                    result_queue,
                    stop_event,
                    log_queue,
                    logger.level,
                    feasible_N,
                    cn_data,
                ),
                name=f'worker-rs-{worker_id + 1}',
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
        stop_requested = run_controller_rs(
            rin,
            cids,
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
        # ---------- cancel task queue join
        task_queue.cancel_join_thread()

        # ---------- terminate workers
        for process in processes:
            if process.is_alive():
                process.terminate()
        for process in processes:
            process.join()

        # ---------- execution error
        if isinstance(error, Exception):
            logger.exception(
                'Random structure generation failed'
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

    # ---------- stop
    if stop_requested:
        raise SystemExit()

    # ---------- time
    if rin.algo == 'RS':
        etime = timedelta(seconds=perf_counter() - time_start)
        logger.info(f'Elapsed time for structure generation: {etime}')
