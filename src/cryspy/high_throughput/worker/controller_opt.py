from enum import Enum, auto
from logging import getLogger
from logging.handlers import QueueListener
import multiprocessing as mp
from pathlib import Path
from queue import Empty
import sqlite3

from ...IO.read_input import ReadInput
from ..db.record import (
    Status,
    claim_next_struc,
    has_record_with_status,
    update_opt_struc,
    update_status,
)
from ..db.sqlite import connect_db
from .worker_opt import run_worker_opt


logger = getLogger('cryspy')


class RunStatus(Enum):
    """High-throughput execution status."""

    COMPLETED = auto()
    STOPPED = auto()


def queue_next_structure(
    conn: sqlite3.Connection,
    task_queue,
) -> int | None:
    """Claim and queue the next structure."""

    # ---------- claim structure
    selected = claim_next_struc(conn)
    if selected is None:
        return None
    record_id, raw_struc_dict = selected

    # ---------- queue task
    task_queue.put(
        (
            record_id,
            raw_struc_dict,
        )
    )

    # ---------- return
    return record_id


def run_controller_opt(
    rin: ReadInput,
    task_queue,
    result_queue,
    stop_event,
    processes: list[mp.Process],
) -> RunStatus:
    """Control database access and worker tasks."""

    # ---------- initialize database connection
    conn = None

    try:
        # ---------- connect database
        conn = connect_db()

        # ---------- initial tasks
        num_active = 0
        no_waiting = False
        stop_requested = False
        worker_failed = False
        for _ in range(rin.njob):
            if Path('STOP_CRYSPY_HT').is_file():
                stop_requested = True
                logger.info(
                    'Stop file detected: STOP_CRYSPY_HT'
                )
                break

            record_id = queue_next_structure(
                conn,
                task_queue,
            )
            if record_id is None:
                no_waiting = True
                break
            num_active += 1

        if num_active == 0:
            if stop_requested:
                return RunStatus.STOPPED

            logger.info('No waiting structures')
            return RunStatus.COMPLETED

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

            record_id, opt_atoms, energy, status = result
            num_active -= 1

            try:
                # ------ register result
                if status == Status.ERROR:
                    update_status(
                        conn,
                        record_id,
                        status,
                    )
                else:
                    update_opt_struc(
                        conn,
                        record_id,
                        opt_atoms,
                        energy,
                        rin.symprec,
                    )
                    update_status(
                        conn,
                        record_id,
                        status,
                    )
                conn.commit()

            except Exception:
                # ------ rollback
                conn.rollback()
                raise

            # ------ finish
            if status != Status.ERROR:
                logger.info(
                    f'    ID {record_id}: '
                    f'E = {energy:.8f} eV/atom, {status.name}'
                )

            # ------ stop assigning tasks
            if worker_failed:
                continue

            # ------ check stop file
            if (
                not stop_requested
                and Path('STOP_CRYSPY_HT').is_file()
            ):
                stop_requested = True
                logger.info(
                    'Stop file detected: STOP_CRYSPY_HT'
                )

            # ------ queue next task
            if not no_waiting and not stop_requested:
                next_id = queue_next_structure(
                    conn,
                    task_queue,
                )
                if next_id is None:
                    no_waiting = True
                else:
                    num_active += 1

        # ---------- worker failure
        if worker_failed:
            raise RuntimeError(
                'Worker stopped before returning all results'
            ) from None

        # ---------- run status
        if stop_requested:
            return RunStatus.STOPPED

        return RunStatus.COMPLETED

    finally:
        # ---------- close database
        if conn is not None:
            conn.close()


def launch_workers_opt(
    rin: ReadInput,
) -> RunStatus:
    """Launch multiple worker processes."""

    # ---------- check stop file
    if Path('STOP_CRYSPY_HT').is_file():
        logger.info(
            'Stop file detected: STOP_CRYSPY_HT'
        )
        return RunStatus.STOPPED

    # ---------- check waiting structures
    conn = connect_db()
    try:
        has_waiting = has_record_with_status(
            conn,
            Status.WAITING,
        )
    finally:
        conn.close()

    if not has_waiting:
        logger.info('No waiting structures')
        return RunStatus.COMPLETED

    # ---------- structure optimization
    logger.info('# ---------- Start structure optimizations')
    logger.info(f'Number of workers: {rin.njob}')
    logger.info(
        f'Minimum interatomic distance check: '
        f'{rin.check_mindist_opt}'
    )

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
        for worker_id in range(rin.njob):
            process = mp.Process(
                target=run_worker_opt,
                args=(
                    rin,
                    task_queue,
                    result_queue,
                    stop_event,
                    log_queue,
                    logger.level,
                ),
                name=f'worker-{worker_id + 1}',
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
        run_status = run_controller_opt(
            rin,
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
                'High-throughput execution failed'
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
                f'{process.name} failed: exitcode = {process.exitcode}'
            )
        raise SystemExit(1)

    # ---------- return
    return run_status
