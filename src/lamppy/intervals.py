"""Candidate target intervals, their uniqueness gate and BED serialisation."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, replace
from pathlib import Path

_EMPTY_FIELD = "."


@dataclass(frozen=True, slots=True)
class Interval:
    """A candidate target region with the metrics the uniqueness gate needs.

    Coordinates are 0-based half-open, matching BED.
    """

    contig: str
    start: int
    end: int
    self_copies_max: int = 1
    genomes_sharing: tuple[str, ...] = ()

    @property
    def name(self) -> str:
        """Stable interval label used as the BED name field."""
        return f"{self.contig}:{self.start}-{self.end}"

    def __len__(self) -> int:
        """Length of the interval in bases."""
        return self.end - self.start


def is_unique(interval: Interval, max_self_copies: int = 1) -> bool:
    """Apply the two-sided uniqueness gate to one interval.

    An interval qualifies only when it is unique intra-genome, occurring at
    most max_self_copies times in its own genome, and unique inter-genome,
    shared with no background genome. Neither condition alone is sufficient:
    a region absent from every background genome but repeated inside its own
    genome yields multi-copy primer binding sites.
    """
    return interval.self_copies_max <= max_self_copies and not interval.genomes_sharing


def filter_unique(
    intervals: Iterable[Interval], max_self_copies: int = 1
) -> list[Interval]:
    """Keep only the intervals that pass the uniqueness gate."""
    return [interval for interval in intervals if is_unique(interval, max_self_copies)]


def drop_short(intervals: Iterable[Interval], min_length: int) -> list[Interval]:
    """Discard intervals shorter than min_length."""
    return [interval for interval in intervals if len(interval) >= min_length]


def merge(intervals: Iterable[Interval], max_gap: int = 0) -> list[Interval]:
    """Merge overlapping or adjacent intervals per contig, combining metrics.

    A merged interval takes the worst intra-genome copy number and the union of
    the sharing genomes, so merging can only make the gate stricter.
    """
    ordered = sorted(intervals, key=lambda item: (item.contig, item.start, item.end))
    merged: list[Interval] = []
    for interval in ordered:
        current = merged[-1] if merged else None
        if (
            current is None
            or current.contig != interval.contig
            or interval.start > current.end + max_gap
        ):
            merged.append(interval)
            continue
        merged[-1] = replace(
            current,
            end=max(current.end, interval.end),
            self_copies_max=max(current.self_copies_max, interval.self_copies_max),
            genomes_sharing=tuple(
                sorted(set(current.genomes_sharing) | set(interval.genomes_sharing))
            ),
        )
    return merged


def write_bed(intervals: Iterable[Interval], path: Path) -> int:
    """Write intervals as BED6 and return how many were written.

    Columns are contig, start, end, name, intra-genome copy number, and the
    comma-separated sharing genomes.
    """
    written = 0
    with path.open("w") as handle:
        for interval in intervals:
            sharing = ",".join(interval.genomes_sharing) or _EMPTY_FIELD
            handle.write(
                f"{interval.contig}\t{interval.start}\t{interval.end}\t"
                f"{interval.name}\t{interval.self_copies_max}\t{sharing}\n"
            )
            written += 1
    return written


def read_bed(path: Path) -> list[Interval]:
    """Read intervals from a BED file written by write_bed."""
    intervals: list[Interval] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            sharing = fields[5] if len(fields) > 5 else _EMPTY_FIELD
            intervals.append(
                Interval(
                    contig=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    self_copies_max=int(fields[4]) if len(fields) > 4 else 1,
                    genomes_sharing=(
                        () if sharing == _EMPTY_FIELD else tuple(sharing.split(","))
                    ),
                )
            )
    return intervals
