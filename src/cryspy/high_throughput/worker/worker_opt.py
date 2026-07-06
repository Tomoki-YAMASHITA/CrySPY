from importlib import util
from logging import getLogger
from logging.handlers import QueueHandler
from math import isfinite
import multiprocessing as mp
from pathlib import Path
from typing import Callable
import warnings

from ase import Atoms

from ...IO.read_input import ReadInput
from ...util.struc_util import check_distance, set_mindist
from ..db.convert import atoms_to_raw, raw_to_atoms, raw_to_struc
from ..db.record import Status


logger = getLogger('cryspy')


class _WarningAggregator:

    _LOGM_PREFIX = (
        'logm result may be inaccurate, approximate err = '
    )

    def __init__(self) -> None:
        self.record_id = None
        self.records = {}

    def start(self, record_id: int) -> None:
        self.record_id = record_id
        self.records.clear()

    def showwarning(
        self,
        message,
        category,
        filename,
        lineno,
        file=None,
        line=None,
    ) -> None:
        text = str(message)

        # ---------- outside structure optimization
        if self.record_id is None:
            logger.warning(
                f'{filename}:{lineno}: '
                f'{category.__name__}: {text}'
            )
            return

        # ---------- aggregate warnings by source
        key = (category.__name__, filename, lineno)
        record = self.records.get(key)
        if record is None:
            record = {
                'count': 0,
                'first_message': text,
                'last_message': text,
                'max_error': None,
            }
            self.records[key] = record

        record['count'] += 1
        record['last_message'] = text

        # ---------- aggregate logm errors
        if text.startswith(self._LOGM_PREFIX):
            try:
                error = float(text[len(self._LOGM_PREFIX):])
            except ValueError:
                error = None

            if error is not None:
                max_error = record['max_error']
                if max_error is None or error > max_error:
                    record['max_error'] = error
            return

        # ---------- log first occurrence
        if record['count'] == 1:
            logger.warning(
                f'ID {self.record_id}: '
                f'{filename}:{lineno}: '
                f'{category.__name__}: {text}'
            )

    def finish(self) -> None:
        for key, record in self.records.items():
            category, filename, lineno = key
            count = record['count']

            # ---------- summarize logm warnings
            if record['first_message'].startswith(
                self._LOGM_PREFIX
            ):
                max_error = record['max_error']
                max_error_text = (
                    'unknown'
                    if max_error is None
                    else f'{max_error:.3e}'
                )
                logger.warning(
                    f'ID {self.record_id}: '
                    f'{filename}:{lineno}: {category}: '
                    f'logm warning occurred {count} times; '
                    f'maximum approximate error = '
                    f'{max_error_text}'
                )

            # ---------- summarize repeated warnings
            elif count > 1:
                logger.warning(
                    f'ID {self.record_id}: '
                    f'{filename}:{lineno}: {category}: '
                    f'warning repeated {count - 1} '
                    f'additional times; last message: '
                    f'{record["last_message"]}'
                )

        self.record_id = None
        self.records.clear()


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

    # ---------- eV/cell --> eV/atom
    energy = energy / len(opt_atoms)

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


def run_worker_opt(
    rin: ReadInput,
    task_queue,
    result_queue,
    stop_event,
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

    # ---------- warning aggregator
    warning_aggregator = _WarningAggregator()
    original_showwarning = warnings.showwarning
    warnings.showwarning = warning_aggregator.showwarning

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
            # ------ stop worker
            if stop_event.is_set():
                break

            # ------ get task
            task = task_queue.get()
            if task is None:
                break

            record_id, raw_struc_dict = task
            logger.info(
                f'{process.name}: '
                f'Start structure optimization: ID {record_id}'
            )
            warning_aggregator.start(record_id)
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

            finally:
                # ------ report aggregated warnings
                warning_aggregator.finish()

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
        # ---------- restore warning handler
        warnings.showwarning = original_showwarning

        # ---------- worker information
        logger.debug(
            f'Finish worker: {process.name}, PID = {process.pid}'
        )
