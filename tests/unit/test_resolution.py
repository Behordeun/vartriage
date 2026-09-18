"""Behaviour of the tri-state lookup-resolution carrier."""

from __future__ import annotations

import pytest

from vartriage.models.resolution import Resolution, ResolutionState


class TestResolutionConstruction:
    def test_found_carries_a_value(self) -> None:
        res = Resolution.found(0.001)
        assert res.state is ResolutionState.FOUND
        assert res.value == 0.001
        assert res.is_found

    def test_confirmed_absent_has_no_value(self) -> None:
        res: Resolution[float] = Resolution.confirmed_absent()
        assert res.state is ResolutionState.CONFIRMED_ABSENT
        assert res.value is None
        assert res.is_confirmed_absent

    def test_lookup_failed_records_a_reason(self) -> None:
        res: Resolution[float] = Resolution.lookup_failed("read timeout")
        assert res.state is ResolutionState.LOOKUP_FAILED
        assert res.value is None
        assert res.reason == "read timeout"
        assert res.is_lookup_failed

    def test_found_rejects_a_none_value(self) -> None:
        with pytest.raises(ValueError):
            Resolution.found(None)


class TestResolutionSemantics:
    def test_confirmed_absent_is_not_a_failure(self) -> None:
        # The whole point of the tri-state: a genuine miss must be
        # distinguishable from a lookup that could not complete.
        absent: Resolution[float] = Resolution.confirmed_absent()
        failed: Resolution[float] = Resolution.lookup_failed("429")
        assert absent.state is not failed.state
        assert not absent.is_lookup_failed
        assert not failed.is_confirmed_absent

    def test_value_or_returns_default_only_when_no_value(self) -> None:
        assert Resolution.found(0.02).value_or(1.0) == 0.02
        assert Resolution.confirmed_absent().value_or(1.0) == 1.0
        assert Resolution.lookup_failed("boom").value_or(1.0) == 1.0

    def test_usable_value_true_only_when_found(self) -> None:
        # A classifier must only key frequency criteria off a usable value;
        # both absent and failed states withhold the value.
        assert Resolution.found(0.0).has_usable_value
        assert not Resolution.confirmed_absent().has_usable_value
        assert not Resolution.lookup_failed("x").has_usable_value

    def test_found_permits_zero_as_a_real_value(self) -> None:
        # A measured 0.0 (e.g. SpliceAI delta of exactly zero) is a real,
        # usable value and must not collapse to absent.
        res = Resolution.found(0.0)
        assert res.is_found
        assert res.value == 0.0
        assert res.has_usable_value

    def test_resolution_is_immutable(self) -> None:
        res = Resolution.found(0.001)
        with pytest.raises((AttributeError, TypeError)):
            res.value = 0.5  # type: ignore[misc]
