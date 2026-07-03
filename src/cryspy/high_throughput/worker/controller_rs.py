from collections.abc import Sequence
from logging import getLogger
from logging.handlers import QueueListener
import multiprocessing as mp
from queue import Empty

from ...IO.read_input import ReadInput
from ..db.record import insert_init_struc
from ..db.sqlite import connect_db
from .worker_rs import run_worker_rs


logger = getLogger('cryspy')


def run_controller_rs(
    rin: ReadInput,
    cids: Sequence[int],
    task_queue,
    result_queue,
    stop_event,
    processes: list[mp.Process],
) -> None:
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
        for _ in processes:
            try:
                cid = next(cid_iter)
            except StopIteration:
                break
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
                # ------ insert structure
                insert_init_struc(
                    conn,
                    cid,
                    struc,
                    rin.atype,
                    rin.symprec,
                )
                conn.commit()

            except Exception:
                # ------ rollback
                conn.rollback()
                raise

            # ------ stop assigning tasks
            if worker_failed:
                continue

            # ------ queue next task
            try:
                next_cid = next(cid_iter)
            except StopIteration:
                continue

            task_queue.put(next_cid)
            num_active += 1

        # ---------- worker failure
        if worker_failed:
            raise RuntimeError(
                'Worker stopped before returning all structures'
            ) from None

    finally:
        # ---------- close database
        if conn is not None:
            conn.close()


def generate_structures_rs(
    rin: ReadInput,
    cids: Sequence[int],
) -> None:
    """Generate and insert random structures."""

    # ---------- nothing to generate
    if not cids:
        return

    # ---------- number of workers
    num_workers = min(rin.njob, len(cids))
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
        run_controller_rs(
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