'''
Data classes for EA offspring generation.
'''

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Callable, Optional, Union

import numpy as np
from pymatgen.core import Structure


@dataclass(frozen=True)
class OffspringTask:
    """Task for generating one offspring."""

    cid: int
    generation: int
    operation: str


@dataclass(frozen=True)
class OffspringResult:
    """Successful offspring generation result."""

    structure: Structure
    parent_ids: tuple[int, ...]
    parent_attempt_count: int


@dataclass(frozen=True)
class TaskFailureResult:
    """Expected offspring generation failure."""

    reason: str
    parent_attempt_count: int


GenerationResult = Union[
    OffspringResult,
    TaskFailureResult,
]

OperationGenerator = Callable[
    [np.random.Generator],
    GenerationResult,
]


def make_task_rng(
    base_seed: int | None,
    cid: int,
    rng: np.random.Generator | None,
) -> np.random.Generator:
    """Make RNG for one offspring task."""

    # ---------- no seed
    if base_seed is None:
        if rng is None:
            return np.random.default_rng()
        return rng

    # ---------- seed
    seed_sequence = np.random.SeedSequence([base_seed, cid])

    # ---------- return
    return np.random.default_rng(seed_sequence)


@dataclass(frozen=True)
class ParentData:
    """Parent structures for offspring generation."""

    structures: Mapping[int, Structure]


@dataclass(frozen=True)
class TournamentSelectionContext:
    """Parameters for tournament selection."""

    parent_ids: tuple[int, ...]
    tournament_size: int


@dataclass(frozen=True)
class RouletteSelectionContext:
    """Parameters for roulette selection."""

    parent_ids: tuple[int, ...]
    cumulative_probabilities: tuple[float, ...]


SelectionContext = Union[
    TournamentSelectionContext,
    RouletteSelectionContext,
]


def build_roulette_selection_context(
    parent_ids: tuple[int, ...],
    fitness: Mapping[int, float],
    a_rlt: float,
    b_rlt: float,
    fit_reverse: bool,
) -> RouletteSelectionContext:
    """Build roulette selection context."""

    # ---------- fitness
    scaled_fitness = np.array(
        [fitness[parent_id] for parent_id in parent_ids],
        dtype=float,
    )

    # ------ reverse
    if not fit_reverse:
        scaled_fitness = -scaled_fitness

    # ------ linear scaling
    a = float(a_rlt)
    b = float(b_rlt)
    fmax = float(scaled_fitness.max())
    fmin = float(scaled_fitness.min())
    if fmax != fmin:
        scaled_fitness = (
            (a - b) / (fmax - fmin) * scaled_fitness
            + (b * fmax - a * fmin) / (fmax - fmin)
        )

    # ---------- cumulative probabilities
    cumulative_probabilities = np.cumsum(
        scaled_fitness / scaled_fitness.sum()
    )

    # ---------- return
    return RouletteSelectionContext(
        parent_ids=parent_ids,
        cumulative_probabilities=tuple(
            float(value)
            for value in cumulative_probabilities
        ),
    )


def select_parent_ids(
    context: SelectionContext,
    n_parent: int,
    rng: np.random.Generator,
) -> tuple[int, ...]:
    """Select parent IDs."""

    # ---------- check input
    if n_parent not in (1, 2):
        raise ValueError('n_parent must be 1 or 2')
    if len(context.parent_ids) < n_parent:
        raise ValueError('Not enough parent candidates')

    # ---------- select parents
    parent_ids = []
    while len(parent_ids) < n_parent:

        # ------ tournament
        if isinstance(context, TournamentSelectionContext):
            tournament_size = min(
                context.tournament_size,
                len(context.parent_ids),
            )
            indices = rng.choice(
                len(context.parent_ids),
                tournament_size,
                replace=False,
            )
            parent_id = context.parent_ids[min(indices)]

        # ------ roulette
        else:
            cumulative = np.asarray(
                context.cumulative_probabilities,
            )
            indices = np.where(cumulative < rng.random())[0]
            selected_index = (
                indices[-1] + 1
                if indices.size != 0
                else 0
            )
            parent_id = context.parent_ids[selected_index]

        # ------ avoid duplicate parents
        if parent_ids and parent_id == parent_ids[0]:
            continue

        parent_ids.append(parent_id)

    # ---------- return
    return tuple(parent_ids)


def generate_with_parent_attempts(
    operation: str,
    selection_context: SelectionContext,
    n_parent: int,
    max_parent_attempts: int,
    generate_attempt: Callable[
        [tuple[int, ...], np.random.Generator],
        Optional[Structure],
    ],
    rng: np.random.Generator,
) -> Union[OffspringResult, TaskFailureResult]:
    """Generate one offspring with parent attempts."""

    # ---------- parent attempts
    for attempt_index in range(max_parent_attempts):
        parent_ids = select_parent_ids(
            selection_context,
            n_parent=n_parent,
            rng=rng,
        )
        offspring = generate_attempt(parent_ids, rng)

        # ------ success
        if offspring is not None:
            return OffspringResult(
                structure=offspring,
                parent_ids=parent_ids,
                parent_attempt_count=attempt_index + 1,
            )

    # ---------- failure
    return TaskFailureResult(
        reason=(
            f'Failed to generate {operation} offspring '
            f'after {max_parent_attempts} parent attempts'
        ),
        parent_attempt_count=max_parent_attempts,
    )


class EAOffspringGenerator:
    """Dispatch offspring generation by operation."""

    def __init__(
        self,
        operation_generators: Mapping[
            str,
            OperationGenerator,
        ],
    ) -> None:
        # ---------- operation generators
        self._operation_generators = dict(
            operation_generators
        )

    def generate_offspring(
        self,
        task: OffspringTask,
        rng: np.random.Generator,
    ) -> GenerationResult:
        """Generate one offspring."""

        # ---------- operation
        try:
            generate = self._operation_generators[
                task.operation
            ]
        except KeyError:
            raise ValueError(
                f'Unsupported operation: {task.operation}'
            ) from None

        # ---------- generate offspring
        return generate(rng)