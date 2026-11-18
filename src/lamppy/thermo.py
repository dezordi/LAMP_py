"""Thermodynamics of single oligos and of a whole LAMP set.

All free energies are kcal/mol, converted from the cal/mol primer3 reports, and
all are negative when a structure is predicted. A LAMP reaction holds four to
six oligos at once, so the set-level numbers are the worst case over every pair
rather than over primer pairs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import combinations, permutations

import primer3

from lamppy.kmers import reverse_complement

_KCAL = 1000.0
_FIVE_PRIME = ("F1c", "B1c")
_OFFSET_PAIRS = (("F2", "F1c"), ("B2", "B1c"))


def gc_fraction(sequence: str) -> float:
    """Return the GC fraction of a sequence, or zero if it is empty."""
    upper = sequence.upper()
    return (upper.count("G") + upper.count("C")) / len(upper) if upper else 0.0


def melting_temperature(sequence: str) -> float:
    """Return the melting temperature of an oligo in degrees Celsius."""
    return float(primer3.calc_tm(sequence.upper()))


def hairpin_dg(sequence: str) -> float:
    """Return the free energy of the most stable hairpin of an oligo."""
    return float(primer3.calc_hairpin(sequence.upper()).dg) / _KCAL


def dimer_dg(first: str, second: str) -> float:
    """Return the free energy of the most stable duplex between two oligos."""
    return float(primer3.calc_heterodimer(first.upper(), second.upper()).dg) / _KCAL


def end_dg(first: str, second: str) -> float:
    """Return the free energy anchored at the 3' end of the first oligo.

    A duplex that leaves a 3' end paired is the one polymerase extends, which is
    why this is scored apart from the most stable duplex overall.
    """
    return float(primer3.calc_end_stability(first.upper(), second.upper()).dg) / _KCAL


def end_stability(sequence: str, five_prime: bool = False, bases: int = 5) -> float:
    """Return the free energy of an oligo's terminal bases against the template.

    End stability is a property of the last few bases rather than of the whole
    duplex: those are what hold the end in place long enough to extend. The 3'
    end is measured by default; F1c and B1c are the regions whose 5' end matters,
    because that is the end that primes once the loop has formed.
    """
    upper = sequence.upper()
    end = upper[:bases] if five_prime else upper[-bases:]
    return dimer_dg(end, reverse_complement(end))


@dataclass(frozen=True, slots=True)
class PrimerThermo:
    """Thermodynamic description of one oligo of a set."""

    name: str
    tm: float
    gc: float
    hairpin: float
    end_stability: float


@dataclass(frozen=True, slots=True)
class SetThermo:
    """Thermodynamic description of a whole LAMP set."""

    primers: tuple[PrimerThermo, ...]
    worst_hairpin: float
    worst_dimer: float
    worst_end_dimer: float
    offsets: Mapping[str, float]

    def primer(self, name: str) -> PrimerThermo | None:
        """Return the description of one oligo, or None if the set lacks it."""
        return next((item for item in self.primers if item.name == name), None)


def evaluate_primer(name: str, sequence: str) -> PrimerThermo:
    """Describe one oligo, measuring the end that matters for its region."""
    return PrimerThermo(
        name=name,
        tm=melting_temperature(sequence),
        gc=gc_fraction(sequence),
        hairpin=hairpin_dg(sequence),
        end_stability=end_stability(sequence, name in _FIVE_PRIME),
    )


def evaluate_set(
    regions: Mapping[str, str], oligos: Mapping[str, str] | None = None
) -> SetThermo:
    """Describe a whole set: its regions one by one, its oligos together.

    Melting temperature and end stability are properties of a region, because
    each region anneals on its own. Hairpins and duplexes are properties of the
    oligos actually pipetted, where F1c and F2 are one molecule, so they are
    measured over those when they are given. Both orientations of every pair are
    measured for 3'-anchored duplexes, because only one of the two ends extends.
    """
    primers = tuple(
        evaluate_primer(name, sequence) for name, sequence in regions.items()
    )
    present = dict(oligos) if oligos else dict(regions)
    pairs = list(combinations(present.values(), 2))
    ordered = list(permutations(present.values(), 2))
    return SetThermo(
        primers=primers,
        worst_hairpin=min(
            (hairpin_dg(sequence) for sequence in present.values()), default=0.0
        ),
        worst_dimer=min((dimer_dg(*pair) for pair in pairs), default=0.0),
        worst_end_dimer=min((end_dg(*pair) for pair in ordered), default=0.0),
        offsets={
            f"{first}_{second}": melting_temperature(regions[second])
            - melting_temperature(regions[first])
            for first, second in _OFFSET_PAIRS
            if first in regions and second in regions
        },
    )
