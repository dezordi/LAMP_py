"""Shared fixtures: synthetic genomes written on demand."""

from __future__ import annotations

import random
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest

GenomeWriter = Callable[[str, Mapping[str, str]], Path]


@pytest.fixture
def write_genome(tmp_path: Path) -> GenomeWriter:
    """Return a factory writing a multi-record FASTA from {contig: sequence}."""

    def factory(name: str, contigs: Mapping[str, str]) -> Path:
        path = tmp_path / f"{name}.fasta"
        with path.open("w") as handle:
            for contig, sequence in contigs.items():
                handle.write(f">{contig} description here\n")
                for start in range(0, len(sequence), 70):
                    handle.write(f"{sequence[start : start + 70]}\n")
        return path

    return factory


@pytest.fixture
def random_sequence() -> Callable[[int, int], str]:
    """Return a factory for reproducible random ACGT sequences."""

    def factory(length: int, seed: int = 0) -> str:
        generator = random.Random(seed)
        return "".join(generator.choice("ACGT") for _ in range(length))

    return factory
