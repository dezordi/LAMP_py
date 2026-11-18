"""Conservation of primer binding sites across the genomes of a target group.

A region shared by every genome of a group is not enough: the assay fails when
one of the six binding sites inside it varies. Presence is tested as an exact
match of the whole oligo on either strand, which is stricter than an alignment
score and needs no multiple alignment.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from lamppy import fasta
from lamppy.kmers import reverse_complement

Group = Mapping[str, tuple[str, ...]]


def load_group(paths: Iterable[Path]) -> dict[str, tuple[str, ...]]:
    """Load the target-group genomes into memory, keyed by file name.

    Sequences are upper-cased, so soft-masking in an input never hides a site.
    """
    return {
        path.name: tuple(sequence.upper() for _, sequence in fasta.read_records(path))
        for path in paths
    }


def holds_site(sequences: Sequence[str], oligo: str) -> bool:
    """Report whether a genome holds an oligo's site on either strand."""
    upper = oligo.upper()
    reverse = reverse_complement(upper)
    return any(upper in sequence or reverse in sequence for sequence in sequences)


def site_conservation(oligo: str, group: Group) -> float:
    """Return the fraction of group genomes that hold an oligo's site.

    An empty group scores one, so a run without a target group is never
    penalised for conservation it was not asked to check.
    """
    if not group:
        return 1.0
    held = sum(1 for sequences in group.values() if holds_site(sequences, oligo))
    return held / len(group)


def conservation_of(oligos: Mapping[str, str], group: Group) -> dict[str, float]:
    """Score every oligo of a set against the target group."""
    return {name: site_conservation(oligo, group) for name, oligo in oligos.items()}
