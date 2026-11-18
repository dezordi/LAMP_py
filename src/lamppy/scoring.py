"""Ranking of candidate primer sets by an explicit, auditable penalty.

A designer returns the first sets it finds, which biases the output toward the
start of the reference. Scoring every candidate and ranking them removes that
bias, and every component of the penalty is reported so a ranking can be
recomputed or overridden without rerunning the pipeline. A high penalty is a
worse set; no candidate is discarded by the score alone.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

from lamppy import primers, profiles, thermo


@dataclass(frozen=True, slots=True)
class Score:
    """A weighted penalty and the components it was summed from."""

    total: float
    components: Mapping[str, float]


@dataclass(frozen=True, slots=True)
class Candidate:
    """A primer set with everything that was measured about it."""

    primer_set: primers.PrimerSet
    measurements: thermo.SetThermo
    conservation: Mapping[str, float]
    score: Score

    @property
    def weakest_site(self) -> float:
        """Conservation of the least conserved binding site of the set."""
        return min(self.conservation.values(), default=1.0)


def _outside(value: float, window: tuple[float, float]) -> float:
    return max(0.0, window[0] - value, value - window[1])


def _below(value: float, limit: float) -> float:
    return max(0.0, limit - value)


def _above(value: float, limit: float) -> float:
    return max(0.0, value - limit)


def components(
    measurements: thermo.SetThermo,
    profile: dict[str, Any],
    conservation: float = 1.0,
) -> dict[str, float]:
    """Measure every penalty component of one set against a profile.

    Each component is zero when the set satisfies the profile and grows with the
    distance from it, in the unit of the parameter it came from.
    """
    values = profiles.parameters(profile)
    offset = values["tm_offset"]["ideal"]
    return {
        "tm": sum(
            _outside(primer.tm, window)
            for primer in measurements.primers
            if (window := profiles.tm_range(profile, primer.name)) is not None
        ),
        "tm_offset": sum(
            abs(value - offset) for value in measurements.offsets.values()
        ),
        "end_stability": sum(
            _above(primer.end_stability, values["end_stability_dg"])
            for primer in measurements.primers
        ),
        "hairpin": _below(measurements.worst_hairpin, values["hairpin_dg"]),
        "dimer": _below(measurements.worst_dimer, values["dimer_dg"]),
        "end_dimer": _below(measurements.worst_end_dimer, values["dimer_dg"]),
        "gc": sum(
            _outside(primer.gc, tuple(values["gc"])) for primer in measurements.primers
        ),
        "conservation": 1.0 - conservation,
    }


def score(
    measurements: thermo.SetThermo,
    profile: dict[str, Any],
    conservation: float = 1.0,
) -> Score:
    """Weight every component of a set's penalty and sum them."""
    parts = components(measurements, profile, conservation)
    return Score(
        total=sum(
            profiles.weight(profile, name) * value for name, value in parts.items()
        ),
        components=parts,
    )


def evaluate(
    primer_set: primers.PrimerSet,
    profile: dict[str, Any],
    conservation: Mapping[str, float] | None = None,
) -> Candidate:
    """Measure and score one primer set."""
    regions = {primer.name: primer.sequence for primer in primer_set.primers}
    measurements = thermo.evaluate_set(regions, primer_set.reaction_oligos())
    sites = dict(conservation) if conservation is not None else {}
    weakest = min(sites.values(), default=1.0)
    return Candidate(
        primer_set=primer_set,
        measurements=measurements,
        conservation=sites,
        score=score(measurements, profile, weakest),
    )


def rank(candidates: Iterable[Candidate]) -> list[Candidate]:
    """Order candidates by increasing penalty, best first."""
    return sorted(
        candidates, key=lambda item: (item.score.total, item.primer_set.index)
    )
