"""Tri-state result for an external-data lookup.

A lookup for a value like an allele frequency, a CADD score, or a gene
constraint metric has three genuinely distinct outcomes, and treating any
two of them as the same produces wrong classifications:

- ``FOUND``: the source was consulted and returned a value.
- ``CONFIRMED_ABSENT``: the source was consulted and the variant/gene is
  not in it. For a population database this is real evidence of rarity.
- ``LOOKUP_FAILED``: the source could not be consulted (timeout, rate
  limit, server error, open circuit, missing/renamed column). This is the
  absence of evidence, not evidence of absence.

Collapsing ``LOOKUP_FAILED`` into ``CONFIRMED_ABSENT`` lets a transient
network failure read downstream as "rare", which biases a variant toward
pathogenic. Collapsing either non-found state into a numeric default has
the same effect. ``Resolution`` keeps the three apart so a consumer can
apply a criterion only on a usable value and withhold it otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Generic, TypeVar

T = TypeVar("T")


class ResolutionState(Enum):
    """Outcome of an external-data lookup."""

    FOUND = "found"
    CONFIRMED_ABSENT = "confirmed_absent"
    LOOKUP_FAILED = "lookup_failed"


@dataclass(frozen=True, slots=True)
class Resolution(Generic[T]):
    """The outcome of one lookup, carrying a value only when ``FOUND``.

    Construct through the ``found``/``confirmed_absent``/``lookup_failed``
    factories rather than the initialiser, so the invariants (a value iff
    found, a reason only on failure) hold by construction.
    """

    state: ResolutionState
    value: T | None = None
    reason: str | None = None

    @classmethod
    def found(cls, value: T) -> Resolution[T]:
        """A consulted source returned ``value``.

        A ``None`` value is rejected: "found nothing" is
        ``confirmed_absent``, not a found ``None``. A falsy-but-real value
        such as ``0.0`` is a valid found value.
        """
        if value is None:
            raise ValueError(
                "Resolution.found requires a value; use confirmed_absent() "
                "for a consulted-but-missing result"
            )
        return cls(state=ResolutionState.FOUND, value=value)

    @classmethod
    def confirmed_absent(cls) -> Resolution[T]:
        """The source was consulted and does not contain the item."""
        return cls(state=ResolutionState.CONFIRMED_ABSENT)

    @classmethod
    def lookup_failed(cls, reason: str) -> Resolution[T]:
        """The source could not be consulted; ``reason`` says why."""
        return cls(state=ResolutionState.LOOKUP_FAILED, reason=reason)

    @property
    def is_found(self) -> bool:
        return self.state is ResolutionState.FOUND

    @property
    def is_confirmed_absent(self) -> bool:
        return self.state is ResolutionState.CONFIRMED_ABSENT

    @property
    def is_lookup_failed(self) -> bool:
        return self.state is ResolutionState.LOOKUP_FAILED

    @property
    def has_usable_value(self) -> bool:
        """True only in the ``FOUND`` state.

        A consumer that applies a criterion from a looked-up value should
        gate on this, so neither a confirmed miss nor a failed lookup is
        mistaken for a real value.
        """
        return self.state is ResolutionState.FOUND

    def value_or(self, default: T) -> T:
        """The found value, or ``default`` when not found.

        The default is returned for both non-found states; a caller that
        needs to treat a confirmed miss differently from a failure must
        branch on ``state`` instead.
        """
        if self.state is ResolutionState.FOUND and self.value is not None:
            return self.value
        return default
