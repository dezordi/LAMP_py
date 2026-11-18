"""The LAMP primer set data model, shared by every designer and consumer.

A set is six named oligos in the coordinates of the sequence they were designed
against, plus the two composites a reaction actually holds. Keeping the model
apart from any one designer is what lets a second backend, the scoring, the
tables and the validation all speak about the same object.
"""

from __future__ import annotations

from dataclasses import dataclass

INNER_PAIRS = (("FIP", "F1c", "F2"), ("BIP", "B1c", "B2"))
INDEPENDENT = ("F3", "B3", "LF", "LB")
ORDER = ("F3", "F2", "F1c", "B1c", "B2", "B3")


@dataclass(frozen=True, slots=True)
class Primer:
    """One oligo of a LAMP set, placed on the genome it was designed against."""

    name: str
    contig: str
    start: int
    sequence: str

    def __len__(self) -> int:
        """Length of the oligo in bases."""
        return len(self.sequence)

    @property
    def end(self) -> int:
        """Position just past the oligo's last base."""
        return self.start + len(self.sequence)


@dataclass(frozen=True, slots=True)
class PrimerSet:
    """A LAMP primer set in the coordinates of its source sequence."""

    reference: str
    index: int
    primers: tuple[Primer, ...]

    def primer(self, name: str) -> Primer | None:
        """Return the oligo with this name, or None if the set lacks it."""
        return next(
            (primer for primer in self.primers if primer.name == name),
            None,
        )

    def composite(self, outer: str, inner: str) -> str | None:
        """Join two oligos into the single one that is ordered, 5' to 3'.

        A LAMP reaction uses FIP, the concatenation of F1c and F2, and BIP, the
        concatenation of B1c and B2, rather than those four oligos separately.
        """
        first, second = self.primer(outer), self.primer(inner)
        if first is None or second is None:
            return None
        return first.sequence + second.sequence

    def reaction_oligos(self) -> dict[str, str]:
        """Return the oligos the reaction actually holds, keyed by name.

        F1c and F2 reach the tube as one molecule, FIP, and B1c and B2 as BIP,
        so interactions between oligos are properties of these and not of the
        six regions a designer reports.
        """
        oligos = {
            name: sequence
            for name, outer, inner in INNER_PAIRS
            if (sequence := self.composite(outer, inner)) is not None
        }
        for name in INDEPENDENT:
            primer = self.primer(name)
            if primer is not None:
                oligos[name] = primer.sequence
        return oligos

    @property
    def span(self) -> tuple[int, int] | None:
        """Outer boundaries of the amplified region, from F3 to B3."""
        first, last = self.primer("F3"), self.primer("B3")
        if first is None or last is None:
            return None
        return first.start, last.end
