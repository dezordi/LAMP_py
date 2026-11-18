"""A LAMP designer in Python, built on primer3-py thermodynamics.

This is the second backend. It exists because no installable LAMP designer
exists besides GLAPD: GLAPD is a Linux-oriented source tarball with no package
in any index, so `pip install lamppy` alone could otherwise never design a set.

The design is deliberately simple and auditable rather than clever. Candidate
sites are enumerated per region with the profile's length and melting
temperature windows, assembled by the geometric constraints of a LAMP set, and
then measured and ranked by the same scoring the GLAPD backend uses. No set is
discarded by the score.

Specificity is not checked here: with this backend it comes from the target
interval itself, which the k-mer gate already proved absent from every
background genome, and from `lamppy validate`, which screens the finished sets
against background, host and database sequence.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from itertools import islice
from typing import Any

from lamppy import fasta, profiles, thermo
from lamppy.kmers import reverse_complement
from lamppy.primers import ORDER, Primer, PrimerSet

_REVERSED = ("F1c", "B2", "B3")
_DEFAULT_LIMIT = 20000
_OVERSAMPLE = 10


@dataclass(frozen=True, slots=True)
class Site:
    """A candidate binding site of one LAMP region on the plus strand."""

    name: str
    start: int
    length: int
    sequence: str
    tm: float

    @property
    def end(self) -> int:
        """Position just past the site's last base on the plus strand."""
        return self.start + self.length


def oligo(sequence: str, start: int, length: int, name: str) -> str:
    """Return the oligo as it is synthesised for one region of a set.

    F1c, B2 and B3 are the regions whose oligo is the reverse complement of the
    plus strand, which is what makes a set's six oligos face each other.
    """
    window = sequence[start : start + length].upper()
    return reverse_complement(window) if name in _REVERSED else window


def sites(sequence: str, name: str, profile: dict[str, Any]) -> list[Site]:
    """Enumerate the candidate sites of one region, best length per position.

    Only the length whose melting temperature is closest to the middle of the
    profile's window is kept per start position: the lengths of one region at
    one position are near-duplicates, and collapsing them is what keeps the
    assembly tractable without hiding a distinct candidate.
    """
    values = profiles.parameters(profile)
    window = profiles.tm_range(profile, name)
    lengths = values["length"].get(name)
    if window is None or lengths is None:
        return []
    ideal = (window[0] + window[1]) / 2
    found: list[Site] = []
    for start in range(len(sequence)):
        candidates = [
            Site(name, start, length, candidate, thermo.melting_temperature(candidate))
            for length in range(lengths[0], lengths[1] + 1)
            if start + length <= len(sequence)
            and (candidate := oligo(sequence, start, length, name)).isalpha()
        ]
        inside = [site for site in candidates if window[0] <= site.tm <= window[1]]
        if inside:
            found.append(min(inside, key=lambda site: abs(site.tm - ideal)))
    return found


def _by_key(found: Sequence[Site], key: str) -> dict[int, list[Site]]:
    index: dict[int, list[Site]] = {}
    for site in found:
        index.setdefault(site.end if key == "end" else site.start, []).append(site)
    return index


def _within(index: Mapping[int, list[Site]], low: int, high: int) -> Iterator[Site]:
    for position in range(low, high + 1):
        yield from index.get(position, ())


def forward_arms(
    sequence: str, profile: dict[str, Any]
) -> Iterator[tuple[Site, Site, Site]]:
    """Yield the best (F3, F2, F1c) triple for each candidate F2 site.

    F3 lies upstream of F2 within the allowed gap, and the 5' end of F1c, which
    is its right edge because its oligo is reverse-complemented, lies the
    required distance downstream of the 5' end of F2.

    Only one triple per F2 site is kept, the one whose melting temperatures sit
    closest to their windows. Keeping every combination is what makes a naive
    enumeration unusable: a 2 kb interval admits hundreds of thousands of them,
    nearly all shifted from one another by a base or two.
    """
    geometry = profiles.parameters(profile)["geometry"]
    outer = _by_key(sites(sequence, "F3", profile), "end")
    inner = _by_key(sites(sequence, "F1c", profile), "end")
    for second in sites(sequence, "F2", profile):
        before = _within(
            outer,
            second.start - geometry["f3_f2_gap"][1],
            second.start - geometry["f3_f2_gap"][0],
        )
        after = [
            site
            for site in _within(
                inner,
                second.start + geometry["f2_f1c"][0],
                second.start + geometry["f2_f1c"][1],
            )
            if site.start >= second.end
        ]
        best = _best_arm(
            ((first, second, last) for first in before for last in after), profile
        )
        if best is not None:
            yield best


def backward_arms(
    sequence: str, profile: dict[str, Any]
) -> Iterator[tuple[Site, Site, Site]]:
    """Yield the best (B1c, B2, B3) triple for each candidate B1c site.

    The constraints mirror the forward ones: the 5' end of B2 is its right edge,
    so the distance from the start of B1c to the end of B2 is the one the
    profile bounds, and B3 follows B2 within the allowed gap.
    """
    geometry = profiles.parameters(profile)["geometry"]
    middle = _by_key(sites(sequence, "B2", profile), "end")
    outer = _by_key(sites(sequence, "B3", profile), "start")
    for first in sites(sequence, "B1c", profile):
        arms = (
            (first, second, third)
            for second in _within(
                middle,
                first.start + geometry["b1c_b2"][0],
                first.start + geometry["b1c_b2"][1],
            )
            if second.start >= first.end
            for third in _within(
                outer,
                second.end + geometry["b2_b3_gap"][0],
                second.end + geometry["b2_b3_gap"][1],
            )
        )
        best = _best_arm(arms, profile)
        if best is not None:
            yield best


def _best_arm(
    arms: Iterator[tuple[Site, Site, Site]] | Any, profile: dict[str, Any]
) -> tuple[Site, Site, Site] | None:
    return min(arms, key=lambda arm: tm_deviation(arm, profile), default=None)


def assemble(
    sequence: str, profile: dict[str, Any], limit: int = _DEFAULT_LIMIT
) -> list[tuple[Site, ...]]:
    """Combine forward and backward arms into complete, geometrically valid sets.

    The number of combinations is bounded, because a long interval admits more
    than any ranking needs; the bound is reported by the caller rather than
    hidden.
    """
    geometry = profiles.parameters(profile)["geometry"]
    gap, amplicon = geometry["f1c_b1c_gap"], geometry["amplicon"]
    backward = _group_by_start(backward_arms(sequence, profile))
    found: list[tuple[Site, ...]] = []
    for first, second, third in forward_arms(sequence, profile):
        for start in range(third.end + gap[0], third.end + gap[1] + 1):
            for fourth, fifth, sixth in backward.get(start, ()):
                span = sixth.end - first.start
                if amplicon[0] <= span <= amplicon[1]:
                    found.append((first, second, third, fourth, fifth, sixth))
                    if len(found) >= limit:
                        return found
    return found


def _group_by_start(
    arms: Iterator[tuple[Site, Site, Site]],
) -> dict[int, list[tuple[Site, Site, Site]]]:
    grouped: dict[int, list[tuple[Site, Site, Site]]] = {}
    for arm in arms:
        grouped.setdefault(arm[0].start, []).append(arm)
    return grouped


def to_primer_set(chosen: Sequence[Site], reference: str, index: int) -> PrimerSet:
    """Turn six sites into a primer set in source genome coordinates."""
    return PrimerSet(
        reference=reference,
        index=index,
        primers=tuple(
            Primer(
                name=name,
                contig=fasta.to_source_coordinate(reference, site.start)[0],
                start=fasta.to_source_coordinate(reference, site.start)[1],
                sequence=site.sequence,
            )
            for name, site in zip(ORDER, chosen, strict=True)
        ),
    )


def tm_deviation(chosen: Sequence[Site], profile: dict[str, Any]) -> float:
    """Total distance of a set's melting temperatures from their ideal windows.

    This is the cheap proxy that decides which assembled sets are worth the full
    thermodynamic scoring: it reuses the temperatures already computed during
    enumeration and costs nothing, while removing the positional bias that
    taking the first sets found would introduce.
    """
    total = 0.0
    for site in chosen:
        window = profiles.tm_range(profile, site.name)
        if window is not None:
            total += abs(site.tm - (window[0] + window[1]) / 2)
    return total


def design_fragment(
    fragment: fasta.Fragment,
    profile: dict[str, Any] | None = None,
    sets: int = 10,
    limit: int = _DEFAULT_LIMIT,
) -> list[PrimerSet]:
    """Design LAMP sets over one target interval and return the candidates.

    More candidates are returned than were asked for, chosen by melting
    temperature alone, so that the full scoring has something to choose between
    and the caller's ranking decides what is reported.
    """
    criteria = profile if profile is not None else profiles.load()
    combinations = assemble(fragment.sequence, criteria, limit)
    combinations.sort(key=lambda chosen: tm_deviation(chosen, criteria))
    return [
        to_primer_set(chosen, fragment.identifier, index)
        for index, chosen in enumerate(
            islice(combinations, sets * _OVERSAMPLE), start=1
        )
    ]
