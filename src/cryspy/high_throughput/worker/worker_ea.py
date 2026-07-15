from logging import getLogger
from logging.handlers import QueueHandler
import multiprocessing as mp
import signal

import numpy as np

from ...EA.ea_offspring import (
    OffspringResult,
    make_task_rng,
)


logger = getLogger('cryspy')


def run_worker_ea(
    offspring_generator,
    base_seed,
    task_queue,
    result_queue,
    stop_event,
    log_queue=None,
    log_level=None,
) -> None:
    """Generate EA offspring without database access."""

    # ---------- ignore SIGINT
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    # ---------- worker logger
    if log_queue is not None:
        logger.handlers.clear()
        logger.propagate = False
        logger.setLevel(log_level)
        logger.addHandler(QueueHandler(log_queue))

    # ---------- worker information
    process = mp.current_process()
    logger.debug(
        f'Start worker: {process.name}, PID = {process.pid}'
    )

    try:
        # ---------- worker RNG for non-seeded runs
        worker_rng = (
            np.random.default_rng()
            if base_seed is None
            else None
        )

        # ---------- generate offspring
        while True:
            # ------ stop worker
            if stop_event.is_set():
                break

            # ------ get task
            task = task_queue.get()
            if task is None:
                break

            logger.info(
                f'{process.name}: '
                f'Start offspring generation: ID {task.cid}, '
                f'operation = {task.operation}'
            )

            try:
                # ------ generate offspring
                rng = make_task_rng(base_seed, task.cid, worker_rng)
                result = offspring_generator.generate_offspring(
                    task,
                    rng,
                )

            except Exception as error:
                # ------ generation error
                logger.exception(
                    f'Offspring generation failed: ID {task.cid}'
                )
                result_queue.put(
                    (
                        task,
                        None,
                        f'{type(error).__name__}: {error}',
                    )
                )
                continue

            # ------ check generation result
            if not isinstance(result, OffspringResult):
                result_queue.put(
                    (
                        task,
                        result,
                        result.reason,
                    )
                )
                continue

            # ------ return result
            result_queue.put(
                (
                    task,
                    result,
                    None,
                )
            )

    except Exception:
        # ---------- worker error
        logger.exception(f'{process.name}: Worker failed')
        raise SystemExit(1)

    finally:
        # ---------- worker information
        logger.debug(
            f'Finish worker: {process.name}, PID = {process.pid}'
        )