from importlib import util
from logging import getLogger
from logging.handlers import QueueHandler
from math import isfinite
import multiprocessing as mp
from pathlib import Path
from typing import Callable

from ase import Atoms

from ...IO.read_input import ReadInput
from ...util.struc_util import check_distance, set_mindist
from ..db.convert import atoms_to_raw, raw_to_atoms, raw_to_struc
from ..db.record import Status


logger = getLogger('cryspy')


def _load_optimize_structure(
    ase_python: str,
) -> Callable[[Atoms], tuple[Atoms, float, bool]]:
    """Load optimize_structure() from calc_in."""

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
    optimize_func = getattr(module, 'optimize_structure', None)
    if not callable(optimize_func):
        raise AttributeError(
            f'optimize_structure() not found in {calc_path}'
        )

    # ---------- return
    return optimize_func


def process_structure(
    rin: ReadInput,
    optimize_func: Callable[[Atoms], tuple[Atoms, float, bool]],
    record_id: int,
    raw_struc_dict: dict,
) -> tuple[Atoms, float, Status]:
    """Process one structure without database access."""

    # ---------- optimize structure
    atoms = raw_to_atoms(raw_struc_dict)
    opt_atoms, energy, converged = optimize_func(atoms)
    if not isinstance(opt_atoms, Atoms):
        raise TypeError(
            'optimize_structure() must return an ASE Atoms object'
        )
    if not isfinite(energy):
        raise ValueError(
            'optimize_structure() must return a finite energy'
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

    # ---------- return
    return opt_atoms, energy, status


def run_worker(
    rin: ReadInput,
    task_queue,
    result_queue,
    log_queue=None,
    log_level=None,
) -> None:
    """Run structure optimizations without database access."""

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
        # ---------- load optimize function
        optimize_func = _load_optimize_structure(rin.ase_python)

        # ---------- run optimization
        while True:
            # ------ get task
            task = task_queue.get()
            if task is None:
                break

            record_id, raw_struc_dict = task
            logger.info(
                f'{process.name}: '
                f'Start structure optimization: ID {record_id}'
            )
            try:
                # ------ optimize structure
                opt_atoms, energy, status = process_structure(
                    rin,
                    optimize_func,
                    record_id,
                    raw_struc_dict,
                )

            except Exception:
                # ------ optimization error
                logger.exception(
                    f'Structure optimization failed: ID {record_id}'
                )
                result_queue.put(
                    (
                        record_id,
                        None,
                        None,
                        Status.ERROR,
                    )
                )
                continue

            # ------ return result
            result_queue.put(
                (
                    record_id,
                    opt_atoms,
                    energy,
                    status,
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
