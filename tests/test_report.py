"""Tests for the HTML report."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from lamppy import report

_SETS = (
    "reference\tset\trank\tcontig\tstart\tend\tamplicon\tscore\n"
    "c1:0\t1\t1\tc1\t100\t320\t220\t4.5\n"
    "c1:0\t2\t2\tc1\t500\t760\t260\t9.1\n"
)
_OLIGOS = (
    "reference\tset\trank\tprimer\tcontig\tstart\tlength\tsequence\n"
    "c1:0\t1\t1\tF3\tc1\t100\t18\tCGGTTCCCTGTACTCGAA\n"
    "c1:0\t1\t1\tF2\tc1\t140\t20\tCGAAATCTTCAAACACTACG\n"
    "c1:0\t1\t1\tF1c\tc1\t180\t23\tTGCTAATTCCTTTGTTGAGGCCG\n"
    "c1:0\t1\t1\tB1c\tc1\t220\t21\tCCTGACGGAAGCAGCATTAAG\n"
    "c1:0\t1\t1\tB2\tc1\t260\t20\tTCGAACGTAACCAAAGTCGT\n"
    "c1:0\t1\t1\tB3\tc1\t302\t18\tACGTTGGCATTCGATTCC\n"
    "c1:0\t1\t1\tFIP\tc1\t180\t43\tTGCTAATTCCTTTGTTGAGGCCGCGAAATCTTCAAACACTACG\n"
    "c1:0\t1\t1\tBIP\tc1\t220\t41\tCCTGACGGAAGCAGCATTAAGTCGAACGTAACCAAAGTCGT\n"
)
_VALIDATION = (
    "set\trole\tsource\tverdict\tsites_found\n"
    "c1:0/1\tbackground\tnear.fasta\tconfiguration\t6\n"
)


@pytest.fixture
def rundir(tmp_path: Path) -> Path:
    directory = tmp_path / "design-a"
    directory.mkdir()
    (directory / "sets.tsv").write_text(_SETS)
    (directory / "primers.tsv").write_text(_OLIGOS)
    (directory / "design.manifest.json").write_text(
        json.dumps(
            {
                "command": "design",
                "lamppy_version": "0.5.0a0",
                "parameters": {
                    "backend": "glapd",
                    "glapd_commit": "0fcce4d",
                    "profile": "eiken-default",
                    "intervals": 1,
                    "num": 10,
                    "weights": {"tm": 1.0, "conservation": 20.0},
                },
                "inputs": [{"path": "/tmp/a.fasta", "sha256": "ab" * 32}],
            }
        )
    )
    return directory


def test_read_table_tolerates_a_missing_file(tmp_path: Path):
    assert report.read_table(tmp_path / "absent.tsv") == []


def test_the_report_names_the_run_and_states_the_caveat(rundir: Path):
    html = report.render(rundir)
    assert "design-a" in html
    assert "wet-lab validation is required" in html
    assert html.startswith("<!doctype html>")


def test_the_report_holds_the_ranked_sets_and_the_manifest(rundir: Path):
    html = report.render(rundir)
    assert "Ranked sets (2)" in html
    assert "ab" * 32 in html
    assert "0.5.0a0" in html


def test_the_order_sheet_lists_the_oligos_that_are_ordered(rundir: Path):
    sheet = report.order_sheet(report.read_table(rundir / "primers.tsv"), "c1:0", "1")
    assert "TGCTAATTCCTTTGTTGAGGCCGCGAAATCTTCAAACACTACG" in sheet
    assert sheet.index("FIP") < sheet.index("BIP")
    assert '<tr><td class="">F1c</td>' not in sheet


def test_the_diagram_is_drawn_to_scale(rundir: Path):
    svg = report.diagram(report.read_table(rundir / "primers.tsv"), "c1:0", "1")
    assert svg.startswith("<svg")
    assert ">100<" in svg
    assert "220 bp" in svg


def test_a_set_without_every_region_is_not_diagrammed(rundir: Path):
    rows = [
        row
        for row in report.read_table(rundir / "primers.tsv")
        if row["primer"] != "B3"
    ]
    assert report.diagram(rows, "c1:0", "1") == ""


def test_only_the_top_sets_are_detailed(rundir: Path):
    assert "Set 2 of c1:0" in report.render(rundir)
    assert "Set 2 of c1:0" not in report.render(rundir, top=1)


def test_validation_evidence_is_included_when_present(rundir: Path):
    (rundir / "validation.tsv").write_text(_VALIDATION)
    html = report.render(rundir)
    assert "Validation (1)" in html
    assert 'class="configuration"' in html


def test_sequences_are_escaped(rundir: Path):
    (rundir / "validation.tsv").write_text(
        "set\tevidence\nc1:0/1\t<script>alert(1)</script>\n"
    )
    assert "<script>" not in report.render(rundir)


def test_write_names_the_file_in_the_run_directory(rundir: Path):
    path = report.write(rundir)
    assert path == rundir / "report.html"
    assert path.read_text().startswith("<!doctype html>")


def test_write_accepts_a_destination(rundir: Path, tmp_path: Path):
    destination = tmp_path / "elsewhere.html"
    assert report.write(rundir, destination) == destination


def test_a_directory_without_tables_is_refused(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="lamppy design"):
        report.write(tmp_path)


def test_column_headings_carry_their_description_as_a_tooltip(rundir: Path):
    html = report.render(rundir)
    assert '<th title="Oligo sequence 5&#x27; to 3&#x27;, as ordered.">sequence' in html


def test_a_column_without_a_description_gets_a_plain_heading():
    assert report.heading("unknown", {}) == "<th>unknown</th>"


def test_every_reported_column_has_a_description():
    from lamppy import tables, validation

    penalties = [
        f"penalty_{name}"
        for name in (
            "tm",
            "tm_offset",
            "end_stability",
            "hairpin",
            "dimer",
            "end_dimer",
            "gc",
            "conservation",
        )
    ]
    fields = [
        *tables._SET_FIELDS,
        *tables._OLIGO_FIELDS[:8],
        *validation._VALIDATION_FIELDS,
        *penalties,
        "dtm_F2_F1c",
        "dtm_B2_B1c",
    ]
    missing = [field for field in fields if field not in report._VALIDATION_COLUMNS]
    assert missing == []


def test_a_set_is_described_by_its_interval_rank_and_amplicon(rundir: Path):
    row = report.read_table(rundir / "sets.tsv")[0]
    text = report.describe_set(row, 2)
    assert "starts at position 0 of contig <code>c1</code>" in text
    assert "Ranked 1 of 2" in text
    assert "c1:100-320, 220 bp" in text


def test_a_reference_without_an_offset_is_named_as_given():
    text = report.describe_set({"reference": "plain", "set": "3"}, 1)
    assert "the target interval <code>plain</code>" in text


def test_the_summary_names_the_designer_its_commit_and_the_genome(rundir: Path):
    html = report.render(rundir)
    assert "<strong>glapd at commit <code>0fcce4d</code></strong>" in html
    assert "<code>a.fasta</code>" in html
    assert "<code>par.pl</code> drops" in html


def test_a_run_without_a_design_manifest_says_the_designer_is_unknown(
    rundir: Path,
):
    (rundir / "design.manifest.json").unlink()
    assert "the designer is unknown" in report.render(rundir)


def test_the_ranking_section_lists_the_weights_the_run_used(rundir: Path):
    html = report.render(rundir)
    assert '<td class="">penalty_conservation</td><td class="">20.0</td>' in html


def test_the_screening_states_its_thresholds_and_tallies_verdicts():
    record = {"parameters": {"min_coverage": 0.9, "max_mismatches": 3, "max_span": 500}}
    text = report.screening(record, [{"verdict": "clear"}, {"verdict": "clear"}])
    assert "at least 90% of the oligo with at most 3 mismatches" in text
    assert "0 configuration, 0 partial, 2 clear" in text


def test_the_order_sheet_names_the_regions_each_oligo_joins(rundir: Path):
    sheet = report.order_sheet(report.read_table(rundir / "primers.tsv"), "c1:0", "1")
    assert '<td class="">FIP</td><td class="">F1c + F2</td>' in sheet


def test_the_diagram_colours_regions_by_the_oligo_that_carries_them(rundir: Path):
    svg = report.diagram(report.read_table(rundir / "primers.tsv"), "c1:0", "1")
    assert svg.count('fill="#2d6cdf"') == 2
    assert "<title>F2, part of FIP</title>" in svg
