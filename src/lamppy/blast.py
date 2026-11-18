"""Short-query BLAST searches and column-wise parsing of its tabular output.

Oligos are too short for default BLAST parameters, so every search here uses
the blastn-short task. Output is read by column position from an explicitly
requested format, never by guessing the default layout.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from lamppy import external, fasta

_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstrand",
    "evalue",
    "bitscore",
)
_FORMAT = f"6 {' '.join(_FIELDS)}"
_PLUS = "plus"


@dataclass(frozen=True, slots=True)
class Hit:
    """One local alignment of an oligo against a subject sequence."""

    query: str
    subject: str
    identity: float
    length: int
    mismatches: int
    gaps: int
    query_start: int
    query_end: int
    start: int
    end: int
    strand: str
    bitscore: float

    @property
    def forward(self) -> bool:
        """Whether the oligo aligned to the plus strand of the subject."""
        return self.strand == _PLUS


def search_command(
    query: Path,
    subject: Path | None = None,
    database: str | None = None,
    threads: int = 1,
    max_targets: int = 5000,
) -> list[str]:
    """Build a blastn-short argument vector against a subject file or a database.

    Raises:
        ValueError: if neither or both of subject and database are given.
    """
    if (subject is None) == (database is None):
        raise ValueError("a search needs exactly one of a subject file or a database")
    target = ["-subject", str(subject)] if subject else ["-db", str(database)]
    return [
        "blastn",
        "-task",
        "blastn-short",
        "-query",
        str(query),
        *target,
        "-outfmt",
        _FORMAT,
        "-max_target_seqs",
        str(max_targets),
        "-num_threads",
        str(threads),
        "-dust",
        "no",
        "-soft_masking",
        "false",
    ]


def parse_hits(output: str) -> Iterator[Hit]:
    """Parse tabular BLAST output into hits with plus-strand coordinates.

    Subject coordinates are normalised so start is always the lower one and
    zero-based, with the strand kept in its own field, which keeps downstream
    ordering tests independent of the alignment direction.
    """
    for line in output.splitlines():
        columns = line.split("\t")
        if len(columns) != len(_FIELDS):
            continue
        subject_start, subject_end = int(columns[8]), int(columns[9])
        yield Hit(
            query=columns[0],
            subject=columns[1],
            identity=float(columns[2]),
            length=int(columns[3]),
            mismatches=int(columns[4]),
            gaps=int(columns[5]),
            query_start=int(columns[6]),
            query_end=int(columns[7]),
            start=min(subject_start, subject_end) - 1,
            end=max(subject_start, subject_end),
            strand=columns[10],
            bitscore=float(columns[12]),
        )


def write_queries(oligos: Mapping[str, str], path: Path) -> int:
    """Write oligos as a FASTA query file keyed by their identifiers."""
    return fasta.write_records(oligos.items(), path)


def search(
    oligos: Mapping[str, str],
    workdir: Path,
    subject: Path | None = None,
    database: str | None = None,
    threads: int = 1,
) -> list[Hit]:
    """Search every oligo against one subject or database and return its hits."""
    workdir.mkdir(parents=True, exist_ok=True)
    query = workdir / "queries.fasta"
    write_queries(oligos, query)
    external.require("blastn")
    output = external.run(
        search_command(query, subject, database, threads=threads),
    )
    return list(parse_hits(output))


def by_subject(hits: Iterable[Hit]) -> dict[str, list[Hit]]:
    """Group hits by subject sequence, each group ordered by position."""
    grouped: dict[str, list[Hit]] = {}
    for hit in hits:
        grouped.setdefault(hit.subject, []).append(hit)
    for group in grouped.values():
        group.sort(key=lambda hit: hit.start)
    return grouped


def best(hits: Sequence[Hit]) -> Hit | None:
    """Return the highest-scoring hit of a group, or None if there is none."""
    return max(hits, key=lambda hit: hit.bitscore, default=None)
