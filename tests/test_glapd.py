"""Tests for the GLAPD stage wiring."""

from __future__ import annotations

from pathlib import Path

import pytest

from lamppy import fasta, glapd


@pytest.fixture
def installation(tmp_path: Path) -> glapd.Glapd:
    root = tmp_path / "GLAPD"
    (root / "Par").mkdir(parents=True)
    for name in ("Single", "LAMP", "par.pl"):
        (root / name).touch()
    return glapd.Glapd.from_root(root)


def test_from_root_rejects_a_directory_without_the_executables(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="not a GLAPD installation"):
        glapd.Glapd.from_root(tmp_path)


def test_from_root_rejects_an_installation_without_parameters(tmp_path: Path):
    root = tmp_path / "GLAPD"
    root.mkdir()
    for name in ("Single", "LAMP", "par.pl"):
        (root / name).touch()
    with pytest.raises(FileNotFoundError, match="Par"):
        glapd.Glapd.from_root(root)


def test_single_command_passes_the_directory_and_parameters(
    installation: glapd.Glapd, tmp_path: Path
):
    command = glapd.single_command(installation, tmp_path / "t.fasta", tmp_path)
    assert command[0] == str(installation.single)
    assert "-dir" in command and str(tmp_path) in command
    assert command[command.index("-par") + 1] == str(installation.parameters)
    assert "-loop" not in command


def test_single_command_adds_the_loop_flag(installation: glapd.Glapd, tmp_path: Path):
    command = glapd.single_command(installation, tmp_path / "t.fasta", tmp_path, True)
    assert "-loop" in command


def test_par_command_runs_perl_with_an_explicit_directory(
    installation: glapd.Glapd, tmp_path: Path
):
    command = glapd.par_command(
        installation,
        tmp_path / "t.fasta",
        tmp_path,
        tmp_path / "index",
        tmp_path / "bg.txt",
        Path("/usr/bin/bowtie"),
    )
    assert command[:2] == ["perl", str(installation.par)]
    assert command[command.index("--dir") + 1] == str(tmp_path)
    assert command[command.index("--specific") + 1] == str(tmp_path / "bg.txt")
    assert command[command.index("--bowtie") + 1] == "/usr/bin/bowtie"
    assert command[command.index("--mis_s") + 1] == "2"


def test_par_command_honours_mismatches_and_threads(
    installation: glapd.Glapd, tmp_path: Path
):
    command = glapd.par_command(
        installation,
        tmp_path / "t.fasta",
        tmp_path,
        tmp_path / "index",
        tmp_path / "bg.txt",
        Path("bowtie"),
        mismatches=3,
        threads=8,
    )
    assert command[command.index("--mis_s") + 1] == "3"
    assert command[command.index("--threads") + 1] == "8"


def test_lamp_command_always_requests_specific_sets(
    installation: glapd.Glapd, tmp_path: Path
):
    command = glapd.lamp_command(installation, tmp_path / "t.fasta", tmp_path, sets=25)
    assert command[0] == str(installation.lamp)
    assert "-specific" in command
    assert command[command.index("-num") + 1] == "25"


def test_contig_names_lists_every_record_of_every_genome(write_genome):
    first = write_genome("g1", {"c1": "ACGT", "c2": "ACGT"})
    second = write_genome("g2", {"c3": "ACGT"})
    assert glapd.contig_names([first, second]) == ["c1", "c2", "c3"]


def test_write_name_list_writes_one_name_per_line(tmp_path: Path):
    path = tmp_path / "bg.txt"
    assert glapd.write_name_list(["c1", "c2"], path) == 2
    assert path.read_text() == "c1\nc2\n"


def test_build_index_rejects_an_empty_background(tmp_path: Path):
    with pytest.raises(ValueError, match="at least one background genome"):
        glapd.build_index([], tmp_path / "index")


def test_workdir_name_is_filesystem_safe():
    fragment = fasta.Fragment("c1", 850, "ACGT")
    assert glapd.workdir_name(fragment) == "c1_850"


def _stage_writer(workdir: Path, calls: list[str], write: bool = True):
    def run(command, cwd=None):
        name = Path(command[0]).name
        calls.append(name)
        if name == "Single":
            for stage in ("Inner", "Outer"):
                (workdir / stage).mkdir(parents=True, exist_ok=True)
                (workdir / stage / "candidate").touch()
        if write and len(command) > 1 and command[1].endswith("par.pl"):
            for stage in ("Inner", "Outer"):
                (workdir / stage / "candidate-specific.txt").touch()
        if name == "LAMP":
            Path(command[command.index("-out") + 1]).touch()
        return ""

    return run


def test_design_fragment_runs_all_three_stages_in_order(
    installation: glapd.Glapd, tmp_path: Path, monkeypatch
):
    calls: list[str] = []
    workdir = tmp_path / "work"
    monkeypatch.setattr(glapd.external, "run", _stage_writer(workdir, calls))
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    fragment = fasta.Fragment("c1", 100, "ACGTACGTAC")
    output = glapd.design_fragment(
        installation, fragment, workdir, tmp_path / "index", ["bg1"]
    )
    assert calls == ["Single", "perl", "LAMP"]
    assert output == workdir / "lamp_primers.txt"
    assert (workdir / "background-list.txt").read_text() == "bg1\n"
    assert (workdir / "target.fasta").read_text() == ">c1:100\nACGTACGTAC\n"


def test_design_fragment_stops_when_par_wrote_no_specificity_files(
    installation: glapd.Glapd, tmp_path: Path, monkeypatch
):
    calls: list[str] = []
    workdir = tmp_path / "work"
    monkeypatch.setattr(
        glapd.external, "run", _stage_writer(workdir, calls, write=False)
    )
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    fragment = fasta.Fragment("c1", 100, "ACGTACGTAC")
    with pytest.raises(RuntimeError, match="no specificity file"):
        glapd.design_fragment(
            installation, fragment, workdir, tmp_path / "index", ["bg1"]
        )
    assert "LAMP" not in calls


def test_require_specificity_files_quotes_the_par_report(tmp_path: Path):
    with pytest.raises(RuntimeError, match="bowtie program"):
        glapd.require_specificity_files(
            tmp_path, "Don't have the bowtie program in /nonexistent/bowtie!"
        )


def test_require_specificity_files_accepts_a_complete_run(tmp_path: Path):
    for stage in ("Inner", "Outer"):
        (tmp_path / stage).mkdir()
        (tmp_path / stage / "candidate-specific.txt").touch()
    glapd.require_specificity_files(tmp_path, "")


_REAL_OUTPUT = """The 1 LAMP primers:
  F3: pos:36,length:18 bp, primer(5'-3'):CGGTTCCCTGTACTCGAA
  F2: pos:74,length:19 bp, primer(5'-3'):AATTCCTTTGTTGAGGCCG
  F1c: pos:115,length:24 bp, primer(5'-3'):CGAAATCTTCAAACACTACGTGCT
  B1c: pos:175,length:22 bp, primer(5'-3'):CCTGACGGAAGCAGCATTAAGT
  B2: pos:237,length:19 bp, primer(5'-3'):CGAACGTAACCAAAGTCGT
  B3: pos:276,length:23 bp, primer(5'-3'):TAAAAAATAAAAAACCGTGCACC
The 2 LAMP primers:
  F3: pos:753,length:23 bp, primer(5'-3'):AGCGAAAAATGAAAATAACGATG
  F2: pos:778,length:25 bp, primer(5'-3'):GTTGGACACGTTTTATTAATTGAAG
  F1c: pos:819,length:25 bp, primer(5'-3'):CGCCAAACCATAATACGTCTTATCA
  B1c: pos:844,length:25 bp, primer(5'-3'):ATTGCCTCTTTATCAGTTCATCCTG
  B2: pos:900,length:20 bp, primer(5'-3'):GACTGTTTTTCTTTCGAAAT
  B3: pos:930,length:23 bp, primer(5'-3'):GGACCACTAAAAATATAAGCATG
"""


@pytest.fixture
def primer_file(tmp_path: Path) -> Path:
    path = tmp_path / "lamp_primers.txt"
    path.write_text(_REAL_OUTPUT)
    return path


def test_parse_primer_sets_reads_every_set(primer_file: Path):
    sets = glapd.parse_primer_sets(primer_file, "chrA:1000")
    assert [primer_set.index for primer_set in sets] == [1, 2]
    assert all(len(primer_set.primers) == 6 for primer_set in sets)


def test_parse_primer_sets_translates_to_source_coordinates(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:1000")[0]
    f3 = first.primer("F3")
    assert f3 is not None
    assert (f3.contig, f3.start, len(f3)) == ("chrA", 1036, 18)


def test_parse_primer_sets_keeps_the_reported_length(primer_file: Path):
    for primer_set in glapd.parse_primer_sets(primer_file, "chrA:0"):
        for primer in primer_set.primers:
            assert len(primer) == len(primer.sequence)


def test_parse_primer_sets_returns_nothing_for_an_empty_file(tmp_path: Path):
    path = tmp_path / "empty.txt"
    path.touch()
    assert glapd.parse_primer_sets(path, "chrA:0") == []


def test_parse_primer_sets_rejects_a_reference_without_an_offset(primer_file: Path):
    with pytest.raises(ValueError, match="fragment identifier"):
        glapd.parse_primer_sets(primer_file, "chrA")


def test_composite_joins_the_inner_primers(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:0")[0]
    assert first.composite("F1c", "F2") == (
        "CGAAATCTTCAAACACTACGTGCT" + "AATTCCTTTGTTGAGGCCG"
    )
    assert first.composite("B1c", "B2") == (
        "CCTGACGGAAGCAGCATTAAGT" + "CGAACGTAACCAAAGTCGT"
    )


def test_composite_returns_none_when_an_oligo_is_absent(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:0")[0]
    assert first.composite("LF", "F2") is None
    assert first.primer("LF") is None


def test_span_covers_f3_to_b3(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:1000")[0]
    assert first.span == (1036, 1299)


def test_reaction_oligos_are_the_molecules_pipetted(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:1000")[0]
    oligos = first.reaction_oligos()
    assert sorted(oligos) == ["B3", "BIP", "F3", "FIP"]
    assert oligos["FIP"] == first.composite("F1c", "F2")
    assert oligos["F3"] == first.primer("F3").sequence


def test_require_candidates_names_the_stage_single_skipped(tmp_path: Path):
    (tmp_path / "Inner").mkdir()
    (tmp_path / "Inner" / "candidate").touch()
    with pytest.raises(RuntimeError, match="no Outer candidate primers"):
        glapd.require_candidates(tmp_path)


def test_design_fragment_stops_when_single_found_no_candidates(
    installation: glapd.Glapd, tmp_path: Path, monkeypatch
):
    calls: list[str] = []
    monkeypatch.setattr(
        glapd.external,
        "run",
        lambda command, cwd=None: calls.append(Path(command[0]).name),
    )
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    fragment = fasta.Fragment("c1", 0, "ACGT")
    with pytest.raises(RuntimeError, match="too short or too low in complexity"):
        glapd.design_fragment(
            installation, fragment, tmp_path / "work", tmp_path / "index", ["bg1"]
        )
    assert calls == ["Single"]


def test_parse_primer_sets_reads_an_empty_file_as_no_sets(tmp_path: Path):
    path = tmp_path / "lamp_primers.txt"
    path.touch()
    assert glapd.parse_primer_sets(path, "c1:0") == []


def test_par_command_adds_the_commonality_list(
    installation: glapd.Glapd, tmp_path: Path
):
    command = glapd.par_command(
        installation,
        tmp_path / "t.fasta",
        tmp_path,
        tmp_path / "index",
        tmp_path / "bg.txt",
        Path("/usr/bin/bowtie"),
        common_list=tmp_path / "group.txt",
        common_mismatches=1,
        loop=True,
    )
    assert command[command.index("--common") + 1] == str(tmp_path / "group.txt")
    assert command[command.index("--mis_c") + 1] == "1"
    assert "--loop" in command


def test_lamp_command_adds_the_common_flag(installation: glapd.Glapd, tmp_path: Path):
    command = glapd.lamp_command(
        installation, tmp_path / "t.fasta", tmp_path, 10, True, True
    )
    assert "-common" in command
    assert "-specific" in command
    assert "-loop" in command


def test_design_fragment_writes_the_commonality_list(
    installation: glapd.Glapd, tmp_path: Path, monkeypatch
):
    commands: list[list[str]] = []

    def run(command, cwd=None):
        commands.append(command)
        name = Path(command[0]).name
        if name == "Single":
            for stage in ("Inner", "Outer"):
                (tmp_path / "work" / stage).mkdir(parents=True, exist_ok=True)
                (tmp_path / "work" / stage / "candidate").touch()
        if len(command) > 1 and command[1].endswith("par.pl"):
            for stage in ("Inner", "Outer"):
                (tmp_path / "work" / stage / "candidate-specific.txt").touch()
        if name == "LAMP":
            Path(command[command.index("-out") + 1]).touch()
        return ""

    monkeypatch.setattr(glapd.external, "run", run)
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    fragment = fasta.Fragment("c1", 0, "ACGT" * 50)
    glapd.design_fragment(
        installation,
        fragment,
        tmp_path / "work",
        tmp_path / "index",
        ["bg1"],
        common=["group1", "group2"],
        common_mismatches=1,
    )
    listing = tmp_path / "work" / "common-list.txt"
    assert listing.read_text().split() == ["group1", "group2"]
    par = next(command for command in commands if command[0] == "perl")
    assert par[par.index("--common") + 1] == str(listing)
    assert "-common" in commands[-1]


def test_span_is_unknown_without_the_outer_primers(primer_file: Path):
    first = glapd.parse_primer_sets(primer_file, "chrA:0")[0]
    partial = glapd.PrimerSet(
        first.reference,
        first.index,
        tuple(primer for primer in first.primers if primer.name != "B3"),
    )
    assert partial.span is None


def test_attempt_returns_the_failure_of_a_stage_that_exits_nonzero(monkeypatch):
    def run(command, cwd=None):
        raise RuntimeError("Single failed with code 152: ")

    monkeypatch.setattr(glapd.external, "run", run)
    assert "code 152" in glapd.attempt(["Single"])


def test_require_primers_accepts_an_empty_file(tmp_path: Path):
    path = tmp_path / "lamp_primers.txt"
    path.touch()
    glapd.require_primers(path, "")


def test_require_primers_fails_when_the_file_is_absent(tmp_path: Path):
    with pytest.raises(RuntimeError, match="wrote no primer file"):
        glapd.require_primers(tmp_path / "absent.txt", "LAMP failed with code 50")


def test_design_fragment_tolerates_a_nonzero_exit_from_every_stage(
    installation: glapd.Glapd, tmp_path: Path, monkeypatch
):
    workdir = tmp_path / "work"

    def run(command, cwd=None):
        name = Path(command[0]).name
        if name == "Single":
            for stage in ("Inner", "Outer"):
                (workdir / stage).mkdir(parents=True, exist_ok=True)
                (workdir / stage / "candidate").touch()
            raise RuntimeError("Single failed with code 152: ")
        if len(command) > 1 and command[1].endswith("par.pl"):
            for stage in ("Inner", "Outer"):
                (workdir / stage / "candidate-specific.txt").touch()
            return ""
        if name == "LAMP":
            Path(command[command.index("-out") + 1]).touch()
            raise RuntimeError("LAMP failed with code 50: ")
        return ""

    monkeypatch.setattr(glapd.external, "run", run)
    monkeypatch.setattr(
        glapd.external, "require", lambda name: Path(f"/usr/bin/{name}")
    )
    output = glapd.design_fragment(
        installation,
        fasta.Fragment("c1", 0, "ACGT" * 50),
        workdir,
        tmp_path / "index",
        ["bg1"],
    )
    assert output == workdir / "lamp_primers.txt"
    assert output.exists()
