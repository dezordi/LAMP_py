"""Invocation of the three GLAPD stages, including the specificity filter.

GLAPD designs single primers with `Single`, filters them against target and
background genomes with `par.pl`, and combines the survivors into LAMP sets
with `LAMP`. Running `Single` and `LAMP` without `par.pl` yields sets that were
never checked against the background, which is why all three stages are wired
here.

One interval is designed per run: `par.pl` concatenates every non-header line
of its reference, so a multi-record reference would place primers at positions
that do not exist in the genome.

The exit status of a GLAPD stage carries no information: `Single` and `LAMP`
return a nonzero code after completing their work normally. Every stage is
therefore judged by the files it should have written, and a failing status is
reported only when those files are absent.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from lamppy import external, fasta
from lamppy.primers import Primer, PrimerSet

logger = logging.getLogger(__name__)

_CANDIDATE = "candidate"
_PRIMERS = "lamp_primers.txt"
_BACKGROUND_LIST = "background-list.txt"
_COMMON_LIST = "common-list.txt"
_REFERENCE = "target.fasta"
_STAGE_DIRECTORIES = ("Inner", "Outer")
_SPECIFIC = "-specific.txt"


@dataclass(frozen=True, slots=True)
class Glapd:
    """Resolved locations of the GLAPD executables and parameter tables."""

    root: Path

    @classmethod
    def from_root(cls, root: Path) -> Glapd:
        """Validate a GLAPD installation directory and return its handle.

        Raises:
            FileNotFoundError: if an executable or the parameter directory is
                missing, which usually means `make` has not been run.
        """
        instance = cls(root.resolve())
        for required in (
            instance.single,
            instance.lamp,
            instance.par,
            instance.parameters,
        ):
            if not required.exists():
                raise FileNotFoundError(f"not a GLAPD installation: {required} missing")
        return instance

    @property
    def commit(self) -> str | None:
        """The source commit recorded beside the build, as the image writes it."""
        marker = self.root / "COMMIT"
        return marker.read_text().strip() if marker.exists() else None

    @property
    def single(self) -> Path:
        """The compiled single-primer finder."""
        return self.root / "Single"

    @property
    def lamp(self) -> Path:
        """The compiled LAMP set builder."""
        return self.root / "LAMP"

    @property
    def par(self) -> Path:
        """The commonality and specificity filter."""
        return self.root / "par.pl"

    @property
    def parameters(self) -> Path:
        """The thermodynamic parameter directory."""
        return self.root / "Par"


def single_command(
    glapd: Glapd, reference: Path, workdir: Path, loop: bool = False
) -> list[str]:
    """Build the argument vector for the Single stage."""
    command = [
        str(glapd.single),
        "-in",
        str(reference),
        "-out",
        _CANDIDATE,
        "-dir",
        str(workdir),
        "-par",
        str(glapd.parameters),
    ]
    if loop:
        command.append("-loop")
    return command


def par_command(
    glapd: Glapd,
    reference: Path,
    workdir: Path,
    index: Path,
    background_list: Path,
    bowtie: Path,
    mismatches: int = 2,
    threads: int = 1,
    loop: bool = False,
    common_list: Path | None = None,
    common_mismatches: int = 0,
) -> list[str]:
    """Build the argument vector for the par.pl specificity stage.

    The directory is always passed explicitly because par.pl defaults to the
    PWD environment variable rather than to its own working directory. A
    commonality list, when given, additionally requires every primer to match
    every genome of the target group.
    """
    command = [
        "perl",
        str(glapd.par),
        "--in",
        _CANDIDATE,
        "--ref",
        str(reference),
        "--dir",
        str(workdir),
        "--specific",
        str(background_list),
        "--bowtie",
        str(bowtie),
        "--index",
        str(index),
        "--mis_s",
        str(mismatches),
        "--threads",
        str(threads),
    ]
    if common_list is not None:
        command += ["--common", str(common_list), "--mis_c", str(common_mismatches)]
    if loop:
        command.append("--loop")
    return command


def lamp_command(
    glapd: Glapd,
    reference: Path,
    workdir: Path,
    sets: int = 10,
    loop: bool = False,
    common: bool = False,
) -> list[str]:
    """Build the argument vector for the LAMP stage.

    The specific flag is always set, so only sets that amplify no background
    genome are reported. The output path is absolute because LAMP applies its
    directory option only to the files it reads, writing the result relative to
    the working directory it inherits.
    """
    command = [
        str(glapd.lamp),
        "-in",
        _CANDIDATE,
        "-ref",
        str(reference),
        "-dir",
        str(workdir),
        "-out",
        str(workdir / _PRIMERS),
        "-num",
        str(sets),
        "-par",
        str(glapd.parameters),
        "-specific",
    ]
    if common:
        command.append("-common")
    if loop:
        command.append("-loop")
    return command


def contig_names(genomes: Iterable[Path]) -> list[str]:
    """List the record identifiers of several FASTA files.

    par.pl matches background genomes by FASTA identifier rather than by file
    name, so every contig of every background genome has to be named.
    """
    return [name for genome in genomes for name, _ in fasta.read_records(genome)]


def write_name_list(names: Iterable[str], path: Path) -> int:
    """Write one record identifier per line and return how many were written."""
    written = 0
    with path.open("w") as handle:
        for name in names:
            handle.write(f"{name}\n")
            written += 1
    return written


def build_index(genomes: Sequence[Path], prefix: Path) -> Path:
    """Build one bowtie index over all background genomes and return its prefix.

    Raises:
        ValueError: if no background genome is given.
    """
    if not genomes:
        raise ValueError("a bowtie index needs at least one background genome")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    external.run(
        [
            str(external.require("bowtie-build")),
            ",".join(str(genome) for genome in genomes),
            str(prefix),
        ]
    )
    return prefix


def design_fragment(
    glapd: Glapd,
    fragment: fasta.Fragment,
    workdir: Path,
    index: Path,
    background: Sequence[str],
    sets: int = 10,
    loop: bool = False,
    mismatches: int = 2,
    threads: int = 1,
    common: Sequence[str] = (),
    common_mismatches: int = 0,
) -> Path:
    """Design LAMP sets for one target interval and return the output file.

    Target-group genomes, when given, are required to carry every primer, so the
    reported sets amplify the whole group and none of the background.

    Raises:
        RuntimeError: if a stage fails, including a par.pl argument error.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    reference = workdir / _REFERENCE
    fasta.write_fragments([fragment], reference)
    listing = workdir / _BACKGROUND_LIST
    write_name_list(background, listing)
    common_list: Path | None = None
    if common:
        common_list = workdir / _COMMON_LIST
        write_name_list(common, common_list)
    bowtie = external.require("bowtie")
    report = attempt(single_command(glapd, reference, workdir, loop))
    require_candidates(workdir, report)
    report = attempt(
        par_command(
            glapd,
            reference,
            workdir,
            index,
            listing,
            bowtie,
            mismatches,
            threads,
            loop,
            common_list,
            common_mismatches,
        )
    )
    require_specificity_files(workdir, report)
    report = attempt(lamp_command(glapd, reference, workdir, sets, loop, bool(common)))
    require_primers(workdir / _PRIMERS, report)
    return workdir / _PRIMERS


def attempt(command: Sequence[str]) -> str:
    """Run a GLAPD stage and return its output, or the failure it reported.

    A nonzero exit status is not an error for these stages, so the message is
    returned for the guard that follows to quote if its files are missing.
    """
    try:
        return external.run(command) or ""
    except RuntimeError as error:
        logger.debug("%s reported a failing status: %s", command[0], error)
        return str(error)


def require_candidates(workdir: Path, report: str = "") -> None:
    """Fail when Single produced no candidate primers for one of the stages.

    Single reports an interval that yields candidates for only one of the two
    stages as a normal run, and the omission would otherwise surface as a
    missing par.pl output, which names the wrong cause.

    Raises:
        RuntimeError: if either candidate file is missing.
    """
    missing = [
        stage
        for stage in _STAGE_DIRECTORIES
        if not (workdir / stage / _CANDIDATE).exists()
    ]
    if missing:
        raise RuntimeError(
            f"Single found no {', '.join(missing)} candidate primers in "
            f"{workdir.name}: the interval is too short or too low in "
            f"complexity{f'. {report.strip()}' if report.strip() else ''}"
        )


def require_specificity_files(workdir: Path, report: str) -> None:
    """Fail when par.pl produced no specificity files.

    par.pl reports argument errors on standard output and still exits zero, so
    the presence of its output files is the only reliable success signal. Left
    unchecked, the run would continue into LAMP with primers that were never
    compared against the background genomes.

    Raises:
        RuntimeError: if either specificity file is missing.
    """
    missing = [
        stage
        for stage in _STAGE_DIRECTORIES
        if not (workdir / stage / f"{_CANDIDATE}{_SPECIFIC}").exists()
    ]
    if missing:
        raise RuntimeError(
            f"par.pl wrote no specificity file for {', '.join(missing)}: "
            f"{report.strip() or 'no output'}"
        )


def require_primers(path: Path, report: str) -> None:
    """Fail when LAMP wrote no output file.

    An empty file is a valid result and means no set survived the filters, but
    no file at all means the stage did not run.

    Raises:
        RuntimeError: if the output file is missing.
    """
    if not path.exists():
        raise RuntimeError(
            f"LAMP wrote no primer file for {path.parent.name}: "
            f"{report.strip() or 'no output'}"
        )


def workdir_name(fragment: fasta.Fragment) -> str:
    """Build a filesystem-safe directory name from a fragment identifier."""
    return fragment.identifier.replace(":", "_")


_SET_HEADER = re.compile(r"^The (\d+) LAMP primers:")
_PRIMER_LINE = re.compile(
    r"^\s*(\S+):\s*pos:(\d+),\s*length:(\d+) bp,\s*primer\(5'-3'\):(\S+)"
)


def parse_primer_sets(path: Path, reference: str) -> list[PrimerSet]:
    """Parse a LAMP output file into primer sets in source genome coordinates.

    GLAPD reports positions relative to the reference it was given, which is one
    target interval, so each position is translated back through the interval
    identifier.

    Raises:
        ValueError: if the reference does not carry a fragment offset.
    """
    sets: list[PrimerSet] = []
    index = 0
    primers: list[Primer] = []
    for line in path.read_text().splitlines():
        header = _SET_HEADER.match(line)
        if header:
            if primers:
                sets.append(PrimerSet(reference, index, tuple(primers)))
            index, primers = int(header.group(1)), []
            continue
        found = _PRIMER_LINE.match(line)
        if found:
            contig, start = fasta.to_source_coordinate(reference, int(found.group(2)))
            primers.append(
                Primer(found.group(1), contig, start, found.group(4).upper())
            )
    if primers:
        sets.append(PrimerSet(reference, index, tuple(primers)))
    return sets
