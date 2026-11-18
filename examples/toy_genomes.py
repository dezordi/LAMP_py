"""Write the three toy genomes the quickstart uses.

The genomes are generated from fixed seeds, so the quickstart is reproducible
without shipping FASTA files in the repository.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

CORE = 40_000
TAIL = 20_000
REPEAT = 5_000
WIDTH = 70


def sequence(length: int, seed: int) -> str:
    """Return a reproducible random ACGT sequence."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def write(path: Path, contig: str, bases: str) -> None:
    """Write one single-contig FASTA."""
    with path.open("w") as handle:
        handle.write(f">{contig}\n")
        for start in range(0, len(bases), WIDTH):
            handle.write(f"{bases[start : start + WIDTH]}\n")


def main(destination: Path) -> None:
    """Write target.fasta, near.fasta and far.fasta into a directory."""
    destination.mkdir(parents=True, exist_ok=True)
    core = sequence(CORE, 1)
    repeat = sequence(REPEAT, 2)
    write(
        destination / "target.fasta",
        "target_c1",
        core + repeat + sequence(TAIL - 2 * REPEAT, 3) + repeat,
    )
    write(destination / "near.fasta", "near_c1", core + sequence(TAIL, 4))
    write(destination / "far.fasta", "far_c1", sequence(CORE + TAIL, 5))
    print(f"wrote 3 genomes -> {destination}")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "genomes"))
