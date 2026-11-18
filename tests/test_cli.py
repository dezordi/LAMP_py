"""Tests for the command-line surface."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from lamppy import __version__, blast, fasta, intervals, manifest, masking
from lamppy.cli import main


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


def test_version_is_reported(runner: CliRunner):
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


@pytest.mark.parametrize(
    "command", ["check-deps", "prepare", "kmers", "targets", "design"]
)
def test_every_subcommand_has_help(runner: CliRunner, command: str):
    result = runner.invoke(main, [command, "--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output


def test_check_deps_lists_every_known_tool(runner: CliRunner):
    from lamppy import external

    result = runner.invoke(main, ["check-deps"])
    assert result.exit_code == 0
    for name in external.TOOLS:
        assert name in result.output


def test_check_deps_reports_glapd_unset(
    runner: CliRunner, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.delenv("GLAPD", raising=False)
    result = runner.invoke(main, ["check-deps"])
    assert "glapd          missing" in result.output


def test_check_deps_rejects_unbuilt_glapd(
    runner: CliRunner, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    monkeypatch.setenv("GLAPD", str(tmp_path))
    result = runner.invoke(main, ["check-deps"])
    assert "not a GLAPD installation" in result.output


def test_prepare_writes_fragments(runner: CliRunner, tmp_path: Path, write_genome):
    genome = write_genome("g", {"c1": "ACGT" * 100 + "N" * 10 + "TTGG" * 100})
    outdir = tmp_path / "out"
    result = runner.invoke(
        main, ["prepare", str(genome), "--outdir", str(outdir), "--min-length", "50"]
    )
    assert result.exit_code == 0
    assert (outdir / "g.fragments.fasta").exists()
    assert "2 fragments" in result.output


def test_prepare_rejects_a_missing_genome(runner: CliRunner, tmp_path: Path):
    result = runner.invoke(
        main, ["prepare", str(tmp_path / "absent.fa"), "--outdir", str(tmp_path)]
    )
    assert result.exit_code != 0


def test_prepare_requires_at_least_one_genome(runner: CliRunner, tmp_path: Path):
    result = runner.invoke(main, ["prepare", "--outdir", str(tmp_path)])
    assert result.exit_code != 0


def test_prepare_reports_a_missing_masking_tool(
    runner: CliRunner, tmp_path: Path, write_genome, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT" * 40})
    monkeypatch.setattr(masking.external, "resolve", lambda name: None)
    result = runner.invoke(
        main,
        ["prepare", str(genome), "--outdir", str(tmp_path / "o"), "--dustmasker"],
    )
    assert result.exit_code != 0
    assert "dustmasker is not on PATH" in result.output


def test_prepare_fragments_the_masked_copy(
    runner: CliRunner, tmp_path: Path, write_genome, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGT" * 100})
    outdir = tmp_path / "o"
    monkeypatch.setattr(masking.external, "require", lambda name: Path(f"/bin/{name}"))

    def run(command, cwd=None):
        if command[0].endswith("dustmasker"):
            outdir.mkdir(parents=True, exist_ok=True)
            fasta.write_records(
                [("c1", "acgt" * 20 + "ACGT" * 80)],
                Path(command[command.index("-out") + 1]),
            )
        return ""

    monkeypatch.setattr(masking.external, "run", run)
    result = runner.invoke(
        main,
        [
            "prepare",
            str(genome),
            "--outdir",
            str(outdir),
            "--dustmasker",
            "--min-length",
            "10",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "masked ->" in result.output
    fragments = dict(fasta.read_records(outdir / "g.fragments.fasta"))
    assert fragments == {"c1:0": "acgt" * 20 + "ACGT" * 80}


def test_kmers_writes_a_sharing_table(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    shared = random_sequence(60, 30)
    first = write_genome("g1", {"c1": shared + random_sequence(60, 31)})
    second = write_genome("g2", {"c1": shared + random_sequence(60, 32)})
    outdir = tmp_path / "out"
    result = runner.invoke(
        main,
        ["kmers", str(first), str(second), "--outdir", str(outdir), "-k", "21"],
    )
    assert result.exit_code == 0
    table = outdir / "kmers.k21.sharing.tsv"
    assert table.exists()
    lines = table.read_text().splitlines()
    assert lines[0].endswith("g1.fasta\tg2.fasta")
    assert f"{len(lines) - 1} shared k-mers" in result.output
    assert (
        json.loads((outdir / "kmers.manifest.json").read_text())["command"] == "kmers"
    )
    assert (outdir / "kmers.k21.spectrum.tsv").exists()
    assert (outdir / "kmers.k21.sharing_matrix.tsv").exists()
    assert (outdir / "kmers.k21.summary.tsv").exists()


def test_targets_writes_one_bed_per_genome(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    first = write_genome("g1", {"c1": random_sequence(300, 33)})
    second = write_genome("g2", {"c1": random_sequence(300, 34)})
    outdir = tmp_path / "out"
    result = runner.invoke(
        main,
        [
            "targets",
            str(first),
            str(second),
            "--outdir",
            str(outdir),
            "-k",
            "21",
            "--min-length",
            "100",
        ],
    )
    assert result.exit_code == 0
    assert (outdir / "g1.targets.bed").read_text() == "c1\t0\t300\tc1:0-300\t1\t.\n"
    assert (outdir / "g2.targets.bed").exists()
    assert "1 targets, 300 bp" in result.output


def test_targets_requires_a_background_genome(
    runner: CliRunner, tmp_path: Path, write_genome
):
    genome = write_genome("g1", {"c1": "ACGT" * 100})
    result = runner.invoke(
        main, ["targets", str(genome), "--outdir", str(tmp_path / "o")]
    )
    assert result.exit_code != 0
    assert "at least two genomes" in result.output


def _design_stages(outdir: Path, calls: list[str]):
    def run(command, cwd=None):
        name = Path(command[0]).name
        calls.append(name)
        if name == "Single":
            workdir = Path(command[command.index("-dir") + 1])
            for stage in ("Inner", "Outer"):
                (workdir / stage).mkdir(parents=True, exist_ok=True)
                (workdir / stage / "candidate").touch()
        if len(command) > 1 and command[1].endswith("par.pl"):
            workdir = Path(command[command.index("--dir") + 1])
            for stage in ("Inner", "Outer"):
                (workdir / stage / "candidate-specific.txt").touch()
        if name == "LAMP":
            Path(command[command.index("-out") + 1]).touch()
        return ""

    return run


def _fake_glapd(root: Path) -> Path:
    (root / "Par").mkdir(parents=True)
    for name in ("Single", "LAMP", "par.pl"):
        (root / name).touch()
    return root


def test_design_runs_every_stage_per_interval(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    from lamppy import glapd

    calls: list[str] = []
    monkeypatch.setattr(
        glapd.external,
        "run",
        _design_stages(tmp_path / "out", calls),
    )
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    genome = write_genome("g1", {"c1": random_sequence(300, 40)})
    background = write_genome("g2", {"c2": random_sequence(300, 41)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(background),
            "--glapd",
            str(_fake_glapd(tmp_path / "GLAPD")),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert calls[:4] == ["bowtie-build", "Single", "perl", "LAMP"]
    assert "c1:0" in result.output


def test_design_rejects_a_bed_without_unique_intervals(
    runner: CliRunner, tmp_path: Path, write_genome
):
    genome = write_genome("g1", {"c1": "ACGT" * 100})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t4\tg2.fasta\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(genome),
            "--glapd",
            str(_fake_glapd(tmp_path / "GLAPD")),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "no unique intervals" in result.output


def test_design_rejects_a_directory_that_is_not_glapd(
    runner: CliRunner, tmp_path: Path, write_genome
):
    genome = write_genome("g1", {"c1": "ACGT" * 100})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(genome),
            "--glapd",
            str(tmp_path),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "not a GLAPD installation" in result.output


def test_design_writes_a_primer_table_and_a_manifest(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    from lamppy import glapd

    outdir = tmp_path / "out"
    calls: list[str] = []
    base = _design_stages(outdir, calls)

    def run(command, cwd=None):
        base(command, cwd)
        if Path(command[0]).name == "LAMP":
            Path(command[command.index("-out") + 1]).write_text(
                "The 1 LAMP primers:\n"
                "  F3: pos:10,length:4 bp, primer(5'-3'):ACGT\n"
                "  F2: pos:30,length:4 bp, primer(5'-3'):TTGG\n"
                "  F1c: pos:50,length:4 bp, primer(5'-3'):CCAA\n"
                "  B1c: pos:70,length:4 bp, primer(5'-3'):GGTT\n"
                "  B2: pos:90,length:4 bp, primer(5'-3'):AATT\n"
                "  B3: pos:110,length:4 bp, primer(5'-3'):CCGG\n"
            )
        return ""

    monkeypatch.setattr(glapd.external, "run", run)
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    genome = write_genome("g1", {"c1": random_sequence(300, 70)})
    background = write_genome("g2", {"c2": random_sequence(300, 71)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(background),
            "--glapd",
            str(_fake_glapd(tmp_path / "GLAPD")),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "1 primer sets" in result.output
    rows = (outdir / "primers.tsv").read_text().splitlines()
    assert len(rows) == 9
    assert rows[1].split("\t")[:6] == ["c1:0", "1", "1", "F3", "c1", "10"]
    sets = (outdir / "sets.tsv").read_text().splitlines()
    header = sets[0].split("\t")
    assert len(sets) == 2
    assert "score" in header
    assert {field for field in header if field.startswith("penalty_")} == {
        "penalty_tm",
        "penalty_tm_offset",
        "penalty_end_stability",
        "penalty_hairpin",
        "penalty_dimer",
        "penalty_end_dimer",
        "penalty_gc",
        "penalty_conservation",
    }
    assert sets[1].split("\t")[header.index("amplicon")] == "104"
    record = json.loads((outdir / "design.manifest.json").read_text())
    assert record["command"] == "design"
    assert record["parameters"]["intervals"] == 1
    assert len(record["inputs"]) == 3
    assert all(len(entry["sha256"]) == 64 for entry in record["inputs"])


_REGIONS = {
    "F3": "CGGTTCCCTGTACTCGAA",
    "F2": "CGAAATCTTCAAACACTACG",
    "F1c": "TGCTAATTCCTTTGTTGAGGCCG",
    "B1c": "CCTGACGGAAGCAGCATTAAG",
    "B2": "TCGAACGTAACCAAAGTCGT",
    "B3": "ACGTTGGCATTCGATTCC",
}


def _primer_table(path: Path) -> Path:
    rows = ["reference\tset\tprimer\tsequence"]
    rows += [f"c1:0\t1\t{name}\t{sequence}" for name, sequence in _REGIONS.items()]
    path.write_text("\n".join(rows) + "\n")
    return path


def test_validate_requires_a_source(runner: CliRunner, tmp_path: Path):
    primers = _primer_table(tmp_path / "primers.tsv")
    result = runner.invoke(
        main, ["validate", "--primers", str(primers), "--outdir", str(tmp_path / "v")]
    )
    assert result.exit_code != 0
    assert "--subject" in result.output


def test_validate_reports_a_background_configuration(
    runner: CliRunner, tmp_path: Path, write_genome, monkeypatch
):
    primers = _primer_table(tmp_path / "primers.tsv")
    background = write_genome("bg", {"c1": "ACGT" * 50})
    strands = {"F3": "plus", "F2": "plus", "F1c": "minus", "B1c": "plus", "B2": "minus"}
    lines = []
    for offset, (name, sequence) in enumerate(_REGIONS.items()):
        start = 100 + offset * 40
        lines.append(
            f"c1:0/1|{name}\tchr9\t100.000\t{len(sequence)}\t0\t0\t1\t{len(sequence)}\t"
            f"{start}\t{start + len(sequence)}\t{strands.get(name, 'minus')}\t"
            f"1e-05\t36.0"
        )

    def run(command, cwd=None):
        return "\n".join(lines) + "\n" if command[0] == "blastn" else ""

    monkeypatch.setattr(blast.external, "run", run)
    monkeypatch.setattr(blast.external, "require", lambda name: Path(f"/bin/{name}"))
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    outdir = tmp_path / "v"
    result = runner.invoke(
        main,
        [
            "validate",
            "--primers",
            str(primers),
            "--subject",
            str(background),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "configuration: 1" in result.output
    rows = (outdir / "validation.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    assert len(rows) == 2
    values = dict(zip(header, rows[1].split("\t"), strict=True))
    assert values["verdict"] == "configuration"
    assert values["role"] == "background"
    assert values["sites_found"] == "6"
    assert values["oriented"] == "True"
    record = json.loads((outdir / "validate.manifest.json").read_text())
    assert record["command"] == "validate"


def test_validate_reports_a_clear_host(
    runner: CliRunner, tmp_path: Path, write_genome, monkeypatch
):
    primers = _primer_table(tmp_path / "primers.tsv")
    host = write_genome("host", {"c1": "ACGT" * 50})
    monkeypatch.setattr(blast.external, "run", lambda command, cwd=None: "")
    monkeypatch.setattr(blast.external, "require", lambda name: Path(f"/bin/{name}"))
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    outdir = tmp_path / "v"
    result = runner.invoke(
        main,
        [
            "validate",
            "--primers",
            str(primers),
            "--host",
            str(host),
            "--outdir",
            str(outdir),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "clear: 1" in result.output


def _lamp_output(path: Path) -> None:
    path.write_text(
        "The 1 LAMP primers:\n"
        "  F3: pos:10,length:18 bp, primer(5'-3'):CGGTTCCCTGTACTCGAA\n"
        "  F2: pos:40,length:20 bp, primer(5'-3'):CGAAATCTTCAAACACTACG\n"
        "  F1c: pos:70,length:23 bp, primer(5'-3'):TGCTAATTCCTTTGTTGAGGCCG\n"
        "  B1c: pos:110,length:21 bp, primer(5'-3'):CCTGACGGAAGCAGCATTAAG\n"
        "  B2: pos:140,length:20 bp, primer(5'-3'):TCGAACGTAACCAAAGTCGT\n"
        "  B3: pos:170,length:18 bp, primer(5'-3'):ACGTTGGCATTCGATTCC\n"
    )


def test_report_refuses_a_directory_without_tables(runner: CliRunner, tmp_path: Path):
    result = runner.invoke(main, ["report", "--run", str(tmp_path)])
    assert result.exit_code != 0
    assert "lamppy design" in result.output


def test_run_executes_every_stage_and_resumes(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    from lamppy import glapd

    calls: list[str] = []
    stages = _design_stages(tmp_path / "out", calls)

    def run(command, cwd=None):
        stages(command, cwd)
        if Path(command[0]).name == "LAMP":
            _lamp_output(Path(command[command.index("-out") + 1]))
        return ""

    for module in (glapd, blast, masking, manifest):
        monkeypatch.setattr(module.external, "run", run)
        monkeypatch.setattr(
            module.external, "require", lambda name: Path(f"/usr/bin/{name}")
        )
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    shared = random_sequence(400, 80)
    target = write_genome("t1", {"c1": shared + random_sequence(600, 81)})
    background = write_genome("b1", {"c2": shared + random_sequence(600, 82)})
    config = tmp_path / "run.yaml"
    glapd_root = _fake_glapd(tmp_path / "GLAPD")
    config.write_text(
        f"targets: [{target}]\nbackground: [{background}]\n"
        f"glapd: {glapd_root}\noutdir: {tmp_path / 'out'}\n"
        "kmer_size: 21\nmin_length: 200\nvalidate: true\n"
    )
    result = runner.invoke(main, ["run", "--config", str(config)])
    assert result.exit_code == 0, result.output
    design_dir = tmp_path / "out" / "design-t1"
    assert (tmp_path / "out" / "targets" / "t1.targets.bed").exists()
    assert (design_dir / "sets.tsv").exists()
    assert (design_dir / "validation.tsv").exists()
    assert (design_dir / "report.html").exists()
    assert "LAMP" in calls
    assert "blastn" in calls

    calls.clear()
    again = runner.invoke(main, ["run", "--config", str(config)])
    assert again.exit_code == 0, again.output
    assert calls == []


def test_validate_reports_a_table_without_the_columns(
    runner: CliRunner, tmp_path: Path
):
    primers = tmp_path / "wrong.tsv"
    primers.write_text("name\tseq\nF3\tACGT\n")
    genome = tmp_path / "bg.fasta"
    genome.write_text(">c1\nACGT\n")
    result = runner.invoke(
        main,
        [
            "validate",
            "--primers",
            str(primers),
            "--subject",
            str(genome),
            "--outdir",
            str(tmp_path / "v"),
        ],
    )
    assert result.exit_code != 0
    assert "'primer' and 'sequence'" in result.output


def test_validate_reports_a_missing_blastn(
    runner: CliRunner, tmp_path: Path, monkeypatch
):
    primers = _primer_table(tmp_path / "primers.tsv")
    genome = tmp_path / "bg.fasta"
    genome.write_text(">c1\nACGT\n")
    monkeypatch.setattr(blast.external, "resolve", lambda name: None)
    result = runner.invoke(
        main,
        [
            "validate",
            "--primers",
            str(primers),
            "--subject",
            str(genome),
            "--outdir",
            str(tmp_path / "v"),
        ],
    )
    assert result.exit_code != 0
    assert "blastn is not on PATH" in result.output


def test_run_reports_a_broken_configuration(runner: CliRunner, tmp_path: Path):
    config = tmp_path / "run.yaml"
    config.write_text("targets: [missing.fasta]\nbackground: [b.fasta]\nglapd: /opt\n")
    result = runner.invoke(main, ["run", "--config", str(config)])
    assert result.exit_code != 0
    assert "missing.fasta" in result.output


def test_design_reports_a_failing_stage_as_an_error(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    from lamppy import glapd

    monkeypatch.setattr(
        glapd.external,
        "run",
        lambda command, cwd=None: (_ for _ in ()).throw(RuntimeError("Single died")),
    )
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    genome = write_genome("g1", {"c1": random_sequence(300, 90)})
    background = write_genome("g2", {"c2": random_sequence(300, 91)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(background),
            "--glapd",
            str(_fake_glapd(tmp_path / "GLAPD")),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "Single died" in result.output
    assert "Traceback" not in result.output


def test_targets_writes_the_region_sequences_for_another_designer(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    shared = random_sequence(400, 60)
    first = write_genome("g1", {"c1": shared + random_sequence(400, 61)})
    second = write_genome("g2", {"c1": shared + random_sequence(400, 62)})
    outdir = tmp_path / "out"
    result = runner.invoke(
        main,
        [
            "targets",
            str(first),
            str(second),
            "--outdir",
            str(outdir),
            "-k",
            "21",
            "--min-length",
            "200",
        ],
    )
    assert result.exit_code == 0, result.output
    records = dict(fasta.read_records(outdir / "g1.targets.fasta"))
    interval = intervals.read_bed(outdir / "g1.targets.bed")[0]
    assert list(records) == [f"c1:{interval.start}"]
    assert len(records[f"c1:{interval.start}"]) == len(interval)


def test_run_without_glapd_stops_after_the_targets(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    shared = random_sequence(400, 63)
    target = write_genome("t1", {"c1": shared + random_sequence(600, 64)})
    background = write_genome("b1", {"c2": shared + random_sequence(600, 65)})
    config = tmp_path / "run.yaml"
    config.write_text(
        f"targets: [{target}]\nbackground: [{background}]\n"
        f"outdir: {tmp_path / 'out'}\nkmer_size: 21\nmin_length: 200\n"
    )
    result = runner.invoke(main, ["run", "--config", str(config)])
    assert result.exit_code == 0, result.output
    assert "no glapd configured" in result.output
    assert (tmp_path / "out" / "targets" / "t1.targets.fasta").exists()
    assert not (tmp_path / "out" / "design-t1").exists()


def test_design_with_the_primer3_backend_needs_no_external_tool(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    genome = write_genome("g1", {"c1": random_sequence(1200, 95)})
    outdir = tmp_path / "out"
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t1200\tc1:0-1200\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(write_genome("g2", {"c2": random_sequence(300, 96)})),
            "--backend",
            "primer3",
            "--outdir",
            str(outdir),
            "--num",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    rows = (outdir / "sets.tsv").read_text().splitlines()
    assert len(rows) == 3
    assert rows[1].split("\t")[2] == "1"
    record = json.loads((outdir / "design.manifest.json").read_text())
    assert record["parameters"]["backend"] == "primer3"
    assert record["parameters"]["glapd"] is None


def test_design_requires_glapd_only_for_the_glapd_backend(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    genome = write_genome("g1", {"c1": random_sequence(300, 97)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(write_genome("g2", {"c2": random_sequence(300, 98)})),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "--glapd is required" in result.output


def test_run_designs_with_the_primer3_backend_without_glapd(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    monkeypatch.setattr(blast.external, "run", lambda command, cwd=None: "")
    monkeypatch.setattr(blast.external, "require", lambda name: Path(f"/bin/{name}"))
    shared = random_sequence(400, 99)
    target = write_genome("t1", {"c1": shared + random_sequence(1200, 100)})
    background = write_genome("b1", {"c2": shared + random_sequence(1200, 101)})
    config = tmp_path / "run.yaml"
    config.write_text(
        f"targets: [{target}]\nbackground: [{background}]\n"
        f"outdir: {tmp_path / 'out'}\nkmer_size: 21\nmin_length: 400\n"
        "backend: primer3\nnum: 2\n"
    )
    result = runner.invoke(main, ["run", "--config", str(config)])
    assert result.exit_code == 0, result.output
    design_dir = tmp_path / "out" / "design-t1"
    assert (design_dir / "sets.tsv").exists()
    assert (design_dir / "validation.tsv").exists()
    assert (design_dir / "report.html").exists()
    assert not (design_dir / "background.1.ebwt").exists()


def test_the_primer3_backend_rejects_loop_primers(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    genome = write_genome("g1", {"c1": random_sequence(1200, 102)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t1200\tc1:0-1200\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--background",
            str(write_genome("g2", {"c2": random_sequence(300, 103)})),
            "--backend",
            "primer3",
            "--loop",
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "does not design loop primers" in result.output


def test_the_glapd_backend_requires_a_background(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence
):
    genome = write_genome("g1", {"c1": random_sequence(300, 104)})
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t300\tc1:0-300\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--glapd",
            str(_fake_glapd(tmp_path / "GLAPD")),
            "--outdir",
            str(tmp_path / "out"),
        ],
    )
    assert result.exit_code != 0
    assert "--background is required" in result.output


def test_the_primer3_backend_designs_without_a_background(
    runner: CliRunner, tmp_path: Path, write_genome, random_sequence, monkeypatch
):
    monkeypatch.setattr(manifest.external, "resolve", lambda name: None)
    genome = write_genome("g1", {"c1": random_sequence(1200, 105)})
    outdir = tmp_path / "out"
    bed = tmp_path / "t.bed"
    bed.write_text("c1\t0\t1200\tc1:0-1200\t1\t.\n")
    result = runner.invoke(
        main,
        [
            "design",
            "--genome",
            str(genome),
            "--targets",
            str(bed),
            "--backend",
            "primer3",
            "--outdir",
            str(outdir),
            "--num",
            "2",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (outdir / "sets.tsv").exists()
