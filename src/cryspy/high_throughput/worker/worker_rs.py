from logging import getLogger
from logging.handlers import QueueHandler
import multiprocessing as mp

from ...IO.read_input import ReadInput
from ...RS.rs_gen import RandomStructureGenerator


logger = getLogger('cryspy')


def run_worker_rs(
    rin: ReadInput,
    log_initialization: bool,
    task_queue,
    result_queue,
    stop_event,
    log_queue=None,
    log_level=None,
) -> None:
    """Generate random structures without database access."""

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
        # ---------- initialize structure generator
        generator = RandomStructureGenerator(
            rin=rin,
            mpi_rank=0,
            log_initialization=log_initialization,
        )

        # ---------- generate structures
        while True:
            # ------ stop worker
            if stop_event.is_set():
                break

            # ------ get task
            cid = task_queue.get()
            if cid is None:
                break

            logger.info(
                f'{process.name}: '
                f'Start structure generation: ID {cid}'
            )

            try:
                # ------ generate structure
                struc = generator.generate(cid)

            except Exception as error:
                # ------ generation error
                logger.exception(
                    f'Structure generation failed: ID {cid}'
                )
                result_queue.put(
                    (
                        cid,
                        None,
                        f'{type(error).__name__}: {error}',
                    )
                )
                continue

            # ------ return result
            result_queue.put(
                (
                    cid,
                    struc,
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