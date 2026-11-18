"""Soft-mask detection and optional repeat-masking of input genomes."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path

from lamppy import external, fasta

logger = logging.getLogger(__name__)

_UPPER = frozenset("ACGT")
_LOWER = frozenset("acgt")
_DEFAULT_SCAN_BYTES = 1 << 20


def detect_soft_masking(path: Path, max_bytes: int = _DEFAULT_SCAN_BYTES) -> bool:
    """Report whether a FASTA file appears soft-masked.

    The signal is mixed case among unambiguous bases within sequence lines. An
    all-lowercase file is a naming convention rather than evidence of masking,
    so it reports False. Ambiguous bases are ignored, because an otherwise
    lowercase assembly with uppercase N-runs is not soft-masked. Headers are
    ignored and the scan stops as soon as both cases are seen.
    """
    seen_upper = False
    seen_lower = False
    scanned = 0
    with path.open() as handle:
        for line in handle:
            scanned += len(line)
            if not line.startswith(">"):
                characters = set(line)
                seen_upper = seen_upper or not characters.isdisjoint(_UPPER)
                seen_lower = seen_lower or not characters.isdisjoint(_LOWER)
                if seen_upper and seen_lower:
                    return True
            if scanned >= max_bytes:
                break
    return False


def warn_if_unmasked(
    paths: Iterable[Path], run_dustmasker: bool, run_trf: bool
) -> list[Path]:
    """Warn once about genomes that are neither soft-masked nor being masked.

    Returns the genomes that triggered the warning. Repetitive regions mostly
    fail the intra-genome uniqueness gate on their own, so this never fails the
    run.
    """
    if run_dustmasker or run_trf:
        return []
    unmasked = [path for path in paths if not detect_soft_masking(path)]
    if unmasked:
        names = ", ".join(path.name for path in unmasked)
        logger.warning(
            "No soft-masking detected in %d genome(s): %s. Repetitive regions may "
            "be reported as candidate targets. Pass --dustmasker and/or --trf to "
            "mask them, or supply soft-masked genomes.",
            len(unmasked),
            names,
        )
    return unmasked


_TRF_PARAMETERS = ("2", "7", "7", "80", "10", "50", "500")
_MASKED_SUFFIX = ".masked.fasta"
_AMBIGUOUS = "N"


def dustmasker_command(tool: Path, genome: Path, destination: Path) -> list[str]:
    """Build the argument vector that soft-masks low-complexity regions."""
    return [
        str(tool),
        "-in",
        str(genome),
        "-outfmt",
        "fasta",
        "-out",
        str(destination),
    ]


def trf_command(tool: Path, genome: Path) -> list[str]:
    """Build the argument vector that masks tandem repeats.

    The match, mismatch and delta weights are trf's documented defaults for
    genomic DNA. HTML output is suppressed and the masked FASTA requested.
    """
    return [str(tool), str(genome), *_TRF_PARAMETERS, "-m", "-h"]


def trf_output(genome: Path, workdir: Path) -> Path:
    """Locate the masked FASTA trf names after its own parameters."""
    suffix = ".".join(_TRF_PARAMETERS)
    return workdir / f"{genome.name}.{suffix}.mask"


def run_dustmasker(genome: Path, workdir: Path) -> Path:
    """Run dustmasker on a genome and return its soft-masked output."""
    destination = workdir / f"{genome.name}.dust.fasta"
    tool = external.require("dustmasker")
    external.run(dustmasker_command(tool, genome, destination))
    return destination


def run_trf(genome: Path, workdir: Path) -> Path:
    """Run trf on a genome and return its repeat-masked output.

    trf writes beside the working directory it is given rather than beside its
    input, and its exit status is not a reliable success signal across
    releases, so the masked file is what the run is judged by.

    Raises:
        RuntimeError: if trf produced no masked file.
    """
    tool = external.require("trf")
    destination = trf_output(genome, workdir)
    try:
        external.run(trf_command(tool, genome), cwd=workdir)
    except RuntimeError:
        if not destination.exists():
            raise
    if not destination.exists():
        raise RuntimeError(f"trf wrote no masked sequence for {genome.name}")
    return destination


def merge_masked(sequence: str, *masked: str) -> str:
    """Lowercase every position a masking tool reported, keeping the original base.

    A tool reports a masked position either in lowercase or as an ambiguous
    base, and the original base is kept so no sequence information is lost.
    Masking from several tools and from the input itself accumulates.

    Raises:
        ValueError: if a masked sequence has a different length.
    """
    for candidate in masked:
        if len(candidate) != len(sequence):
            raise ValueError(
                f"masked sequence of {len(candidate)} bases does not match the "
                f"{len(sequence)} bases of the input"
            )
    return "".join(
        base.lower() if any(_is_masked(other[index]) for other in masked) else base
        for index, base in enumerate(sequence)
    )


def _is_masked(base: str) -> bool:
    return base.islower() or base.upper() == _AMBIGUOUS


def soft_mask(
    genome: Path, outdir: Path, run_dust: bool = False, run_repeats: bool = False
) -> Path:
    """Write a soft-masked copy of a genome and return its path.

    The input is never modified. With no tool requested the copy is a faithful
    one, so callers can treat the returned path as the genome to work from
    either way.

    Raises:
        ValueError: if a tool dropped or renamed a record of the input.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    if run_dust:
        outputs.append(run_dustmasker(genome, outdir))
    if run_repeats:
        outputs.append(run_trf(genome, outdir))
    masked = [dict(fasta.read_records(path)) for path in outputs]
    destination = outdir / f"{genome.stem}{_MASKED_SUFFIX}"
    fasta.write_records(_merged_records(genome, outputs, masked), destination)
    return destination


def _merged_records(
    genome: Path, outputs: Sequence[Path], masked: Sequence[Mapping[str, str]]
) -> Iterator[tuple[str, str]]:
    for name, sequence in fasta.read_records(genome):
        for path, record in zip(outputs, masked, strict=True):
            if name not in record:
                raise ValueError(f"{path.name} holds no record named {name!r}")
        yield name, merge_masked(sequence, *(record[name] for record in masked))
