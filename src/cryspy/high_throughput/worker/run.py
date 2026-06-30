from importlib import util
from logging import getLogger
from logging.handlers import QueueHandler, QueueListener
from math import isfinite
import multiprocessing as mp
from pathlib import Path
import sqlite3
from typing import Callable, Optional

from ase import Atoms

from ...IO.read_input import ReadInput
from ...util.struc_util import check_distance, set_mindist
from ..db.convert import atoms_to_raw, raw_to_atoms, raw_to_struc
from ..db.record import (
    Status,
    claim_next_struc,
    update_opt_struc,
    update_status,
)
from ..db.sqlite import connect_db


logger = getLogger('cryspy')


def _load_optimize_atoms(
    ase_python: str,
) -> Callable[[Atoms], tuple[Atoms, float, bool]]:
    """Load optimize_atoms() from calc_in."""

    # ---------- file path
    calc_path = Path('calc_in') / ase_python

    # ---------- load module
    spec = util.spec_from_file_location(
        'cryspy_ht_user_calculator',
        calc_path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f'Could not load {calc_path}')
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # ---------- optimize function
    optimize_atoms = getattr(module, 'optimize_atoms', None)
    if not callable(optimize_atoms):
        raise AttributeError(
            f'optimize_atoms() not found in {calc_path}'
        )

    # ---------- return
    return optimize_atoms


def run_one(
    rin: ReadInput,
    optimize_atoms: Callable[[Atoms], tuple[Atoms, float, bool]],
    conn: sqlite3.Connection,
) -> Optional[int]:
    """Optimize and register one structure."""

    # ---------- claim structure
    selected = claim_next_struc(conn)
    if selected is None:
        logger.info(
            f'{mp.current_process().name}: No waiting structures'
        )
        return None
    record_id, raw_struc_dict = selected
    logger.debug(
        f'{mp.current_process().name} claimed ID {record_id}'
    )
    logger.info(f'Start structure optimization: ID {record_id}')

    try:
        # ---------- optimize structure
        atoms = raw_to_atoms(raw_struc_dict)
        opt_atoms, energy, converged = optimize_atoms(atoms)
        if not isinstance(opt_atoms, Atoms):
            raise TypeError(
                'optimize_atoms() must return an ASE Atoms object'
            )
        if not isfinite(energy):
            raise ValueError(
                'optimize_atoms() must return a finite energy'
            )

        # ---------- check mindist
        mindist_ok = True
        if rin.check_mindist_opt:
            opt_struc = raw_to_struc(atoms_to_raw(opt_atoms))
            mindist = set_mindist(
                rin.atype,
                rin.mindist,
                rin.mindist_factor,
                rin.struc_mode,
                no_print=True,
            )
            mindist_ok, mindist_ij, dist = check_distance(
                opt_struc,
                rin.atype,
                mindist,
            )
            if not mindist_ok:
                type0 = rin.atype[mindist_ij[0]]
                type1 = rin.atype[mindist_ij[1]]
                logger.warning(
                    f'ID {record_id}: '
                    f'mindist: {type0} - {type1}, {dist}'
                )

        # ---------- status
        if not mindist_ok:
            status = Status.MINDIST
        elif not converged:
            status = Status.NOT_CONVERGED
        else:
            status = Status.DONE

        # ---------- register result
        update_opt_struc(
            conn,
            record_id,
            opt_atoms,
            energy,
        )
        update_status(conn, record_id, status)
        conn.commit()

    except Exception:
        # ---------- rollback
        conn.rollback()
        # ---------- update status
        update_status(conn, record_id, Status.ERROR)
        conn.commit()
        logger.exception(
            f'Structure optimization failed: ID {record_id}'
        )
        # ---------- return
        return record_id

    # ---------- finish
    logger.info(
        f'    ID {record_id}: '
        f'E = {energy:.8f} eV/cell, {status.name}'
    )

    # ---------- return
    return record_id


def run_worker(
    rin: ReadInput,
    log_queue=None,
    log_level=None,
) -> None:
    """Run optimizations until no waiting structures remain."""

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

    # ---------- connect database
    conn = connect_db()

    try:
        # ---------- load optimize function
        optimize_atoms = _load_optimize_atoms(rin.ase_python)

        # ---------- run optimization
        while True:
            # ------ check stop file
            if Path('STOP_CRYSPY_HT').is_file():
                logger.info(
                    f'{process.name}: Stop file detected: STOP_CRYSPY_HT'
                )
                break

            # ------ run one optimization
            record_id = run_one(rin, optimize_atoms, conn)
            if record_id is None:
                break

    finally:
        # ---------- close database
        conn.close()

        # ---------- worker information
        logger.debug(
            f'Finish worker: {process.name}, PID = {process.pid}'
        )


def launch_workers(
    rin: ReadInput,
) -> None:
    """Launch multiple worker processes."""

    # ---------- structure optimization
    logger.info('# ---------- Start structure optimizations')
    logger.info(f'Number of workers: {rin.njob}')
    logger.info(
        f'Minimum interatomic distance check: '
        f'{rin.check_mindist_opt}'
    )

    # ---------- single worker
    if rin.njob == 1:
        run_worker(rin)
        return

    # ---------- logging
    log_queue = mp.Queue()
    listener = QueueListener(
        log_queue,
        *logger.handlers,
        respect_handler_level=True,
    )

    # ---------- start workers
    processes = []
    for worker_id in range(rin.njob):
        process = mp.Process(
            target=run_worker,
            args=(rin, log_queue, logger.level),
            name=f'worker-{worker_id + 1}',
        )
        process.start()
        processes.append(process)
        logger.debug(
            f'Launched {process.name}: PID = {process.pid}'
        )

    # ---------- start listener
    listener.start()

    try:
        # ---------- wait workers
        for process in processes:
            process.join()
    finally:
        # ---------- stop logging
        listener.stop()
        log_queue.close()

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
