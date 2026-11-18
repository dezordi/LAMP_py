"""Derivation of candidate target intervals from k-mer occurrence counts."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from itertools import groupby
from pathlib import Path

from lamppy import fasta, kmers, tables, thermo
from lamppy.intervals import Interval, drop_short


def _is_soft_masked(window: str) -> bool:
    return window != window.upper()


def qualifying_positions(
    sequence: str,
    size: int,
    self_counts: Mapping[str, int],
    background: Mapping[str, Counter[str]],
    max_self_copies: int = 1,
) -> dict[int, int]:
    """Return {offset: intra-genome copy number} for windows passing the gate.

    A window qualifies only when its canonical k-mer occurs at most
    max_self_copies times in its own genome and in no background genome at all.
    Soft-masked windows are skipped, since lowercase marks repeat content the
    caller has already flagged.
    """
    qualifying: dict[int, int] = {}
    for offset, kmer in kmers.iter_kmers(sequence, size):
        copies = self_counts.get(kmer, 0)
        if copies > max_self_copies or _is_soft_masked(
            sequence[offset : offset + size]
        ):
            continue
        if any(kmer in counts for counts in background.values()):
            continue
        qualifying[offset] = copies
    return qualifying


def _runs(offsets: list[int]) -> list[list[int]]:
    return [
        [offset for _, offset in group]
        for _, group in groupby(enumerate(offsets), lambda pair: pair[1] - pair[0])
    ]


def fragment_intervals(
    fragment: fasta.Fragment,
    size: int,
    self_counts: Mapping[str, int],
    background: Mapping[str, Counter[str]],
    max_self_copies: int = 1,
) -> list[Interval]:
    """Collapse the qualifying windows of a fragment into source-coordinate intervals.

    Consecutive qualifying windows become one interval spanning them, carrying
    the worst intra-genome copy number observed inside it.
    """
    qualifying = qualifying_positions(
        fragment.sequence, size, self_counts, background, max_self_copies
    )
    return [
        Interval(
            contig=fragment.contig,
            start=fragment.start + run[0],
            end=fragment.start + run[-1] + size,
            self_copies_max=max(qualifying[offset] for offset in run),
        )
        for run in _runs(sorted(qualifying))
    ]


def find_targets(
    genome: Path,
    counts: Mapping[str, Counter[str]],
    size: int,
    min_length: int = 200,
    max_self_copies: int = 1,
) -> list[Interval]:
    """Find the intervals of one genome that are unique against all the others.

    Args:
        genome: The target genome, which must be a key of counts by file name.
        counts: Canonical k-mer counts per genome, including the target itself.
        size: The k-mer length the counts were built with.
        min_length: Discard intervals shorter than this many bases.
        max_self_copies: Tolerated intra-genome copy number.

    Raises:
        KeyError: if counts holds no entry for the target genome.
    """
    self_counts = counts[genome.name]
    background = {
        label: genome_counts
        for label, genome_counts in counts.items()
        if label != genome.name
    }
    found = [
        interval
        for fragment in fasta.iter_fragments(genome, min_length=size)
        for interval in fragment_intervals(
            fragment, size, self_counts, background, max_self_copies
        )
    ]
    return drop_short(found, min_length)


def extract(genome: Path, selected: Iterable[Interval]) -> list[fasta.Fragment]:
    """Cut the sequence of each interval out of its genome.

    The fragment keeps the interval start as its offset, so a position reported
    by a primer designer translates back to the source contig.

    Raises:
        KeyError: if an interval names a contig absent from the genome.
    """
    sequences = dict(fasta.read_records(genome))
    return [
        fasta.Fragment(
            contig=interval.contig,
            start=interval.start,
            sequence=sequences[interval.contig][interval.start : interval.end],
        )
        for interval in selected
    ]


_METRIC_FIELDS = (
    "contig",
    "start",
    "end",
    "length",
    "gc",
    "self_copies_max",
    "genomes_sharing",
)


def metrics(genome: Path, selected: Iterable[Interval]) -> list[dict[str, object]]:
    """Describe each interval so a human can review it before ordering primers.

    GC content is the metric that decides feasibility before any design runs: a
    region too AT-rich for the melting temperatures a LAMP set needs will yield
    no sets however unique it is.
    """
    intervals = list(selected)
    return [
        {
            "contig": fragment.contig,
            "start": interval.start,
            "end": interval.end,
            "length": len(interval),
            "gc": thermo.gc_fraction(fragment.sequence),
            "self_copies_max": interval.self_copies_max,
            "genomes_sharing": ",".join(interval.genomes_sharing) or ".",
        }
        for interval, fragment in zip(
            intervals, extract(genome, intervals), strict=True
        )
    ]


def write_metrics(genome: Path, selected: Iterable[Interval], path: Path) -> int:
    """Write the per-interval metrics table and return how many rows it holds."""
    return tables.write_rows(path, _METRIC_FIELDS, metrics(genome, selected))
