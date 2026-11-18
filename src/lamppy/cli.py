"""Command-line entry points."""

from __future__ import annotations

import logging
import os
from pathlib import Path

import click

from lamppy import (
    __version__,
    blast,
    conservation,
    external,
    fasta,
    glapd,
    intervals,
    kmers,
    manifest,
    masking,
    native,
    pipeline,
    profiles,
    report,
    scoring,
    tables,
    targets,
    validation,
)
from lamppy import (
    primers as primer_model,
)

_GLAPD = "glapd"
_PRIMER3 = "primer3"

_GENOMES = click.argument(
    "genomes",
    nargs=-1,
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
_OUTDIR = click.option(
    "--outdir",
    required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="Directory for the output files.",
)
_KMER_SIZE = click.option(
    "-k",
    "--kmer-size",
    default=31,
    show_default=True,
    help="Canonical k-mer length.",
)


def _configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(levelname)s %(message)s",
    )


@click.group()
@click.version_option(__version__, prog_name="lamppy")
@click.option("-v", "--verbose", is_flag=True, help="Emit debug-level logging.")
def main(verbose: bool) -> None:
    """Design group-specific LAMP primer sets from whole genomes."""
    _configure_logging(verbose)


@main.command("check-deps")
def check_deps() -> None:
    """Report which external tools are installed and their versions."""
    for name in external.TOOLS:
        path = external.resolve(name)
        if path is None:
            click.echo(f"{name:<14} missing")
            continue
        version = external.tool_version(name)
        click.echo(f"{name:<14} {path} {version or ''}".rstrip())
    root = os.environ.get("GLAPD")
    if root is None:
        click.echo(f"{'glapd':<14} missing (set GLAPD or pass --glapd)")
        return
    try:
        glapd.Glapd.from_root(Path(root))
    except FileNotFoundError as error:
        click.echo(f"{'glapd':<14} {error}")
        return
    click.echo(f"{'glapd':<14} {root}")


@main.command()
@_GENOMES
@_OUTDIR
@click.option(
    "--min-length",
    default=200,
    show_default=True,
    help="Discard fragments shorter than this many bases.",
)
@click.option("--dustmasker", is_flag=True, help="Mask low-complexity regions.")
@click.option("--trf", is_flag=True, help="Mask tandem repeats.")
def prepare(
    genomes: tuple[Path, ...],
    outdir: Path,
    min_length: int,
    dustmasker: bool,
    trf: bool,
) -> None:
    """Split genomes into unambiguous fragments, leaving the inputs untouched.

    Ambiguity runs are never spanned, so no fragment crosses an assembly gap.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    masking.warn_if_unmasked(genomes, dustmasker, trf)
    for genome in genomes:
        source = genome
        if dustmasker or trf:
            try:
                source = masking.soft_mask(genome, outdir, dustmasker, trf)
            except (FileNotFoundError, RuntimeError, ValueError) as error:
                raise click.ClickException(str(error)) from error
            click.echo(f"{genome.name}: masked -> {source}")
        destination = outdir / f"{genome.stem}.fragments.fasta"
        count = fasta.write_fragments(
            fasta.iter_fragments(source, min_length), destination
        )
        click.echo(f"{genome.name}: {count} fragments -> {destination}")
    manifest.write(
        outdir,
        "prepare",
        genomes,
        {"min_length": min_length, "dustmasker": dustmasker, "trf": trf},
    )


@main.command("kmers")
@_GENOMES
@_OUTDIR
@_KMER_SIZE
@click.option(
    "--min-genomes",
    default=2,
    show_default=True,
    help="Report k-mers seen in at least this many genomes.",
)
def kmers_command(
    genomes: tuple[Path, ...], outdir: Path, kmer_size: int, min_genomes: int
) -> None:
    """Count canonical k-mers and describe how they are shared across genomes.

    Writes a spectrum, a pairwise sharing matrix and a per-genome summary, which
    together say whether a genome set has private sequence to target at all, plus
    the full table of shared k-mers with one occurrence-count column per genome.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    counts = kmers.count_genomes(genomes, kmer_size)
    prefix = outdir / f"kmers.k{kmer_size}"
    writers = (
        ("spectrum", kmers.write_spectrum),
        ("sharing_matrix", kmers.write_sharing_matrix),
        ("summary", kmers.write_summary),
    )
    for name, writer in writers:
        destination = prefix.with_name(f"{prefix.name}.{name}.tsv")
        click.echo(f"{writer(counts, destination)} rows -> {destination}")
    destination = prefix.with_name(f"{prefix.name}.sharing.tsv")
    rows = kmers.write_sharing_table(counts, destination, min_genomes)
    click.echo(f"{rows} shared k-mers -> {destination}")
    for genomes_held, held in kmers.spectrum(counts).items():
        click.echo(f"k-mers in {genomes_held} genome(s): {held}")
    manifest.write(
        outdir, "kmers", genomes, {"kmer_size": kmer_size, "min_genomes": min_genomes}
    )


@main.command("targets")
@_GENOMES
@_OUTDIR
@_KMER_SIZE
@click.option(
    "--min-length",
    default=200,
    show_default=True,
    help="Discard target intervals shorter than this many bases.",
)
@click.option(
    "--max-self-copies",
    default=1,
    show_default=True,
    help="Tolerated intra-genome copy number of a target k-mer.",
)
def targets_command(
    genomes: tuple[Path, ...],
    outdir: Path,
    kmer_size: int,
    min_length: int,
    max_self_copies: int,
) -> None:
    """Find per-genome regions that are unique intra-genome and inter-genome.

    Every other input genome is the background of the genome being processed,
    so a region survives only when it is single-copy in its own genome and
    absent from all the others. The regions are written as BED, as a metrics
    table and as FASTA, so a run can stop here and feed any primer designer.
    """
    if len(genomes) < 2:
        raise click.ClickException("at least two genomes are needed for a background")
    outdir.mkdir(parents=True, exist_ok=True)
    counts = kmers.count_genomes(genomes, kmer_size)
    for genome in genomes:
        found = targets.find_targets(
            genome, counts, kmer_size, min_length, max_self_copies
        )
        destination = outdir / f"{genome.stem}.targets.bed"
        written = intervals.write_bed(found, destination)
        targets.write_metrics(
            genome, found, outdir / f"{genome.stem}.targets.metrics.tsv"
        )
        sequences = outdir / f"{genome.stem}.targets.fasta"
        fasta.write_fragments(targets.extract(genome, found), sequences)
        total = sum(len(interval) for interval in found)
        click.echo(f"{genome.name}: {written} targets, {total} bp -> {destination}")
        click.echo(f"{genome.name}: sequences -> {sequences}")
    manifest.write(
        outdir,
        "targets",
        genomes,
        {
            "kmer_size": kmer_size,
            "min_length": min_length,
            "max_self_copies": max_self_copies,
        },
    )


@main.command("design")
@click.option(
    "--genome",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="The target genome the primers must amplify.",
)
@click.option(
    "--targets",
    "targets_bed",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="BED file of unique target intervals, as written by 'lamppy targets'.",
)
@click.option(
    "--background",
    "background",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A genome the primers must not amplify; repeatable, and required by the "
    "glapd backend, which aligns against it.",
)
@click.option(
    "--backend",
    type=click.Choice([_GLAPD, _PRIMER3]),
    default=_GLAPD,
    show_default=True,
    help="Which designer to use. The primer3 backend needs no external tool.",
)
@click.option(
    "--glapd",
    "glapd_root",
    envvar="GLAPD",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="The GLAPD installation directory, after running make; "
    "required by the glapd backend.",
)
@click.option(
    "--group",
    "group",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A further genome of the target group the primers must also amplify; "
    "repeatable.",
)
@_OUTDIR
@click.option(
    "--profile",
    default=profiles.DEFAULT,
    show_default=True,
    help="Design-criteria profile: a bundled name or a path to a YAML file.",
)
@click.option(
    "--num", default=10, show_default=True, help="LAMP sets to report per interval."
)
@click.option(
    "--loop", is_flag=True, help="Also design loop primers; glapd backend only."
)
@click.option(
    "--mismatches",
    default=2,
    show_default=True,
    help="Mismatches tolerated when matching a primer to a background genome.",
)
@click.option("--threads", default=1, show_default=True, help="Threads for bowtie.")
@click.option(
    "--mis-c",
    "common_mismatches",
    default=0,
    show_default=True,
    help="Mismatches tolerated when matching a primer to a target-group genome.",
)
@click.option(
    "--max-self-copies",
    default=1,
    show_default=True,
    help="Tolerated intra-genome copy number of a target interval.",
)
def design(
    genome: Path,
    targets_bed: Path,
    background: tuple[Path, ...],
    backend: str,
    glapd_root: Path | None,
    group: tuple[Path, ...],
    outdir: Path,
    profile: str,
    num: int,
    loop: bool,
    mismatches: int,
    threads: int,
    common_mismatches: int,
    max_self_copies: int,
) -> None:
    """Design LAMP primer sets for unique target intervals.

    With the glapd backend all three GLAPD stages run, so every reported set has
    also been checked against the background genomes by par.pl. The primer3
    backend designs in Python and needs no external tool, which makes its sets
    dependent on 'lamppy validate' for a background check. Either way every set
    is measured, scored against the profile and ranked, so the output is not the
    first sets the designer happened to find.
    """
    if backend == _GLAPD and glapd_root is None:
        raise click.ClickException("--glapd is required by the glapd backend")
    if backend == _GLAPD and not background:
        raise click.ClickException("--background is required by the glapd backend")
    if backend == _PRIMER3 and loop:
        raise click.ClickException(
            "the primer3 backend does not design loop primers yet; "
            "use --backend glapd for --loop"
        )
    try:
        criteria = profiles.load(profile)
        installation = (
            glapd.Glapd.from_root(glapd_root) if glapd_root is not None else None
        )
    except (FileNotFoundError, ValueError) as error:
        raise click.ClickException(str(error)) from error
    selected = intervals.filter_unique(intervals.read_bed(targets_bed), max_self_copies)
    if not selected:
        raise click.ClickException(f"no unique intervals in {targets_bed}")
    outdir.mkdir(parents=True, exist_ok=True)
    target_group = conservation.load_group([genome, *group])
    designed: list[primer_model.PrimerSet] = []
    try:
        if installation is not None:
            index = glapd.build_index([*background, *group], outdir / "background")
            names = glapd.contig_names(background)
            common = glapd.contig_names(group)
        for fragment in targets.extract(genome, selected):
            if installation is not None:
                primers = glapd.design_fragment(
                    installation,
                    fragment,
                    outdir / glapd.workdir_name(fragment),
                    index,
                    names,
                    num,
                    loop,
                    mismatches,
                    threads,
                    common,
                    common_mismatches,
                )
                sets = glapd.parse_primer_sets(primers, fragment.identifier)
            else:
                sets = native.design_fragment(fragment, criteria, num)
            designed.extend(sets)
            click.echo(f"{fragment.identifier}: {len(sets)} primer sets")
    except (FileNotFoundError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error
    ranked = scoring.rank(
        scoring.evaluate(
            primer_set,
            criteria,
            conservation.conservation_of(
                {primer.name: primer.sequence for primer in primer_set.primers},
                target_group,
            ),
        )
        for primer_set in designed
    )[: num * len(selected)]
    oligos, sets_table = outdir / "primers.tsv", outdir / "sets.tsv"
    click.echo(f"{tables.write_oligos(ranked, oligos)} oligos -> {oligos}")
    click.echo(f"{tables.write_sets(ranked, sets_table)} ranked sets -> {sets_table}")
    manifest.write(
        outdir,
        "design",
        [genome, targets_bed, *background, *group],
        {
            "backend": backend,
            "glapd": str(installation.root) if installation else None,
            "glapd_commit": installation.commit if installation else None,
            "profile": criteria.get("name", profile),
            "weights": profiles.parameters(criteria)["weights"],
            "intervals": len(selected),
            "num": num,
            "loop": loop,
            "mismatches": mismatches,
            "common_mismatches": common_mismatches,
            "group": [path.name for path in group],
            "max_self_copies": max_self_copies,
        },
    )


@main.command("validate")
@click.option(
    "--primers",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="TSV of oligos, as written by 'lamppy design' or typed by hand.",
)
@click.option(
    "--subject",
    "subjects",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A background genome the sets must not amplify; repeatable.",
)
@click.option(
    "--host",
    "hosts",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="A host or matrix genome to screen against; repeatable.",
)
@click.option(
    "--db",
    "databases",
    multiple=True,
    help="A BLAST database to screen against; repeatable.",
)
@_OUTDIR
@click.option(
    "--max-span",
    default=500,
    show_default=True,
    help="Largest distance from F3 to B3 that still counts as amplifiable.",
)
@click.option(
    "--max-mismatches",
    default=3,
    show_default=True,
    help="Mismatches an alignment may carry and still count as a binding site.",
)
@click.option(
    "--min-coverage",
    default=0.9,
    show_default=True,
    help="Fraction of an oligo an alignment must cover to count as a binding site.",
)
@click.option("--threads", default=1, show_default=True, help="Threads for blastn.")
def validate_command(
    primers: Path,
    subjects: tuple[Path, ...],
    hosts: tuple[Path, ...],
    databases: tuple[str, ...],
    outdir: Path,
    max_span: int,
    max_mismatches: int,
    min_coverage: float,
    threads: int,
) -> None:
    """Screen finished LAMP sets against background, host and database sequence.

    Each oligo is searched with blastn-short, and every source is then tested
    for the arrangement that amplification requires: all six sites in order and
    within one amplicon. The report is evidence per set and source, not a
    verdict of specific, and never replaces wet-lab validation.
    """
    sources: list[tuple[str, str, Path | None]] = [
        *(("background", path.name, path) for path in subjects),
        *(("host", path.name, path) for path in hosts),
        *(("database", name, None) for name in databases),
    ]
    if not sources:
        raise click.ClickException("give at least one --subject, --host or --db")
    try:
        sets = validation.read_primer_table(primers)
    except ValueError as error:
        raise click.ClickException(str(error)) from error
    outdir.mkdir(parents=True, exist_ok=True)
    oligos = {
        name: sequence
        for primer_set in sets
        for name, sequence in validation.queries(primer_set).items()
    }
    rows: list[dict[str, object]] = []
    for role, name, path in sources:
        try:
            hits = blast.search(
                oligos,
                outdir / "blast",
                subject=path,
                database=None if path else name,
                threads=threads,
            )
        except (FileNotFoundError, RuntimeError) as error:
            raise click.ClickException(str(error)) from error
        rows.extend(
            validation.evaluate(
                primer_set, hits, role, name, max_span, max_mismatches, min_coverage
            )
            for primer_set in sets
        )
        click.echo(f"{name}: {len(hits)} hits over {len(sets)} set(s)")
    report = outdir / "validation.tsv"
    click.echo(f"{validation.write_report(rows, report)} rows -> {report}")
    for state in (validation.CONFIGURATION, validation.PARTIAL, validation.CLEAR):
        click.echo(f"{state}: {sum(1 for row in rows if row['verdict'] == state)}")
    manifest.write(
        outdir,
        "validate",
        [primers, *subjects, *hosts],
        {
            "sets": len(sets),
            "databases": list(databases),
            "max_span": max_span,
            "max_mismatches": max_mismatches,
            "min_coverage": min_coverage,
        },
    )


@main.command("report")
@click.option(
    "--run",
    "rundir",
    required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="A run directory written by 'lamppy design' or 'lamppy validate'.",
)
@click.option(
    "--out",
    "destination",
    type=click.Path(dir_okay=False, path_type=Path),
    help="Where to write the HTML; defaults to report.html in the run directory.",
)
@click.option(
    "--top",
    default=5,
    show_default=True,
    help="Sets to detail with a diagram and an order sheet.",
)
def report_command(rundir: Path, destination: Path | None, top: int) -> None:
    """Render one HTML report over the tables of a run directory."""
    try:
        path = report.write(rundir, destination, top)
    except FileNotFoundError as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"report -> {path}")


@main.command("run")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="YAML configuration of a whole pipeline run.",
)
@click.option("--force", is_flag=True, help="Rerun stages whose output already exists.")
@click.pass_context
def run_command(context: click.Context, config_path: Path, force: bool) -> None:
    """Run targets, design, validation and reporting from one configuration.

    Every stage is the command of the same name, so a configured run and a
    hand-driven one take the same code path. A stage whose output is already in
    place is skipped unless --force is given, which makes a run resumable
    without a workflow engine. With the glapd backend and no glapd path the run
    stops after the target regions, which need no external tool at all.
    """
    try:
        settings = pipeline.load(config_path)
    except (ValueError, FileNotFoundError) as error:
        raise click.ClickException(str(error)) from error
    outdir = Path(settings["outdir"])
    for genome in settings["targets"]:
        stem = Path(genome).stem
        beds, design_dir = outdir / "targets", outdir / f"design-{stem}"
        bed = beds / f"{stem}.targets.bed"
        if force or not bed.exists():
            context.invoke(
                targets_command,
                genomes=(genome, *settings["background"]),
                outdir=beds,
                kmer_size=settings["kmer_size"],
                min_length=settings["min_length"],
                max_self_copies=settings["max_self_copies"],
            )
        if settings["backend"] == _GLAPD and settings["glapd"] is None:
            click.echo(
                f"{stem}: no glapd configured, stopping after the target regions"
            )
            continue
        if force or not (design_dir / "sets.tsv").exists():
            context.invoke(
                design,
                genome=genome,
                targets_bed=bed,
                background=settings["background"],
                backend=settings["backend"],
                glapd_root=settings["glapd"],
                group=settings["group"],
                outdir=design_dir,
                profile=settings["profile"],
                num=settings["num"],
                loop=settings["loop"],
                mismatches=settings["mismatches"],
                threads=settings["threads"],
                common_mismatches=settings["common_mismatches"],
                max_self_copies=settings["max_self_copies"],
            )
        if settings["validate"] and (
            force or not (design_dir / "validation.tsv").exists()
        ):
            context.invoke(
                validate_command,
                primers=design_dir / "primers.tsv",
                subjects=settings["background"],
                hosts=settings["host"],
                databases=settings["db"],
                outdir=design_dir,
                max_span=settings["max_span"],
                threads=settings["threads"],
            )
        context.invoke(report_command, rundir=design_dir, destination=None, top=5)
