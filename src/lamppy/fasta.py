"""Read-only FASTA access and splitting of contigs into unambiguous fragments."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

_UNAMBIGUOUS = re.compile(r"[ACGTacgt]+")
_ID_SEPARATOR = ":"


@dataclass(frozen=True, slots=True)
class Fragment:
    """An unambiguous slice of a contig, carrying its offset in that contig.

    The offset is kept in the identifier so positions reported by downstream
    tools translate back to the source contig without a side-car mapping.
    """

    contig: str
    start: int
    sequence: str

    @property
    def identifier(self) -> str:
        """FASTA identifier encoding the source contig and the start offset."""
        return f"{self.contig}{_ID_SEPARATOR}{self.start}"

    def __len__(self) -> int:
        """Length of the fragment in bases."""
        return len(self.sequence)


def read_records(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (identifier, sequence) per record, preserving case.

    Case is preserved because it carries soft-masking information. The
    identifier is truncated at the first whitespace, matching how BLAST,
    bowtie and GenMap report sequence names.
    """
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def split_at_ambiguous(
    contig: str, sequence: str, min_length: int = 1
) -> list[Fragment]:
    """Split a sequence on runs of ambiguous bases into unambiguous fragments.

    Any character outside ACGT (either case) is a break, so N-runs and IUPAC
    ambiguity codes are never spanned. Joining across such a run would create a
    junction absent from the genome, and a primer set designed over it cannot
    amplify. Fragments shorter than min_length are discarded.
    """
    return [
        Fragment(contig, match.start(), match.group())
        for match in _UNAMBIGUOUS.finditer(sequence)
        if len(match.group()) >= min_length
    ]


def iter_fragments(path: Path, min_length: int = 1) -> Iterator[Fragment]:
    """Yield unambiguous fragments for every record in a FASTA file."""
    for contig, sequence in read_records(path):
        yield from split_at_ambiguous(contig, sequence, min_length)


def write_records(
    records: Iterable[tuple[str, str]], path: Path, width: int = 60
) -> int:
    """Write (identifier, sequence) pairs as wrapped FASTA and return the count."""
    written = 0
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(f"{sequence[start : start + width]}\n")
            written += 1
    return written


def write_fragments(fragments: Iterable[Fragment], path: Path, width: int = 60) -> int:
    """Write fragments as FASTA and return how many were written."""
    return write_records(
        ((fragment.identifier, fragment.sequence) for fragment in fragments),
        path,
        width,
    )


def to_source_coordinate(identifier: str, position: int) -> tuple[str, int]:
    """Translate a position within a fragment back to (contig, contig position).

    Raises:
        ValueError: if the identifier does not carry a fragment offset.
    """
    contig, separator, offset = identifier.rpartition(_ID_SEPARATOR)
    if not separator or not offset.isdigit():
        raise ValueError(f"not a fragment identifier: {identifier!r}")
    return contig, int(offset) + position
