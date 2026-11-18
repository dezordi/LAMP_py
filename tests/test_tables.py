"""Tests for the final oligo and set tables."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from lamppy import glapd, profiles, scoring, tables

PRIMERS = {
    "F3": "CGGTTCCCTGTACTCGAA",
    "F2": "CGAAATCTTCAAACACTACG",
    "F1c": "TGCTAATTCCTTTGTTGAGGCCG",
    "B1c": "CCTGACGGAAGCAGCATTAAG",
    "B2": "TCGAACGTAACCAAAGTCGT",
    "B3": "ACGTTGGCATTCGATTCC",
}
_ORDER = ("F3", "F2", "F1c", "B1c", "B2", "B3")


def _primer_set(index: int = 1) -> glapd.PrimerSet:
    return glapd.PrimerSet(
        reference="c1:1000",
        index=index,
        primers=tuple(
            glapd.Primer(name, "c1", 1000 + offset * 40, PRIMERS[name])
            for offset, name in enumerate(_ORDER)
        ),
    )


@pytest.fixture
def candidates() -> list[scoring.Candidate]:
    profile = profiles.load()
    return [
        scoring.evaluate(
            _primer_set(index), profile, dict.fromkeys(PRIMERS, 1.0 - 0.1 * index)
        )
        for index in (1, 2)
    ]


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_the_oligo_table_holds_the_regions_and_the_composites(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "primers.tsv"
    assert tables.write_oligos(candidates, path) == 16
    rows = _read(path)
    names = [row["primer"] for row in rows[:8]]
    assert names == [*_ORDER, "FIP", "BIP"]


def test_an_oligo_row_carries_its_own_thermodynamics(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "primers.tsv"
    tables.write_oligos(candidates, path)
    row = next(row for row in _read(path) if row["primer"] == "F1c")
    assert row["contig"] == "c1"
    assert row["start"] == "1080"
    assert float(row["tm"]) > 50.0
    assert float(row["end_stability_dg"]) < 0.0
    assert float(row["conservation"]) == 0.9


def test_a_composite_row_carries_the_length_that_is_ordered(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "primers.tsv"
    tables.write_oligos(candidates, path)
    row = next(row for row in _read(path) if row["primer"] == "FIP")
    assert row["sequence"] == PRIMERS["F1c"] + PRIMERS["F2"]
    assert int(row["length"]) == len(PRIMERS["F1c"]) + len(PRIMERS["F2"])
    assert row["start"] == "1080"
    assert row["tm"] == ""


def test_the_oligo_table_reports_the_rank_of_its_set(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "primers.tsv"
    tables.write_oligos(scoring.rank(candidates), path)
    ranks = {row["set"]: row["rank"] for row in _read(path)}
    assert ranks == {"1": "1", "2": "2"}


def test_the_set_table_holds_one_row_per_set_with_its_span(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "sets.tsv"
    assert tables.write_sets(candidates, path) == 2
    row = _read(path)[0]
    assert row["reference"] == "c1:1000"
    assert row["start"] == "1000"
    assert row["end"] == str(1200 + len(PRIMERS["B3"]))
    assert int(row["amplicon"]) == 200 + len(PRIMERS["B3"])


def test_every_penalty_component_is_a_column(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "sets.tsv"
    tables.write_sets(candidates, path)
    row = _read(path)[0]
    components = candidates[0].score.components
    for name, value in components.items():
        assert float(row[f"penalty_{name}"]) == pytest.approx(value, abs=1e-3)
    assert float(row["score"]) == pytest.approx(candidates[0].score.total, abs=1e-3)


def test_the_set_table_reports_the_offsets_and_the_worst_interactions(
    candidates: list[scoring.Candidate], tmp_path: Path
):
    path = tmp_path / "sets.tsv"
    tables.write_sets(candidates, path)
    row = _read(path)[0]
    assert float(row["dtm_F2_F1c"]) > 0.0
    assert float(row["worst_dimer_dg"]) < 0.0
    assert float(row["worst_end_dimer_dg"]) < 0.0
    assert float(row["min_conservation"]) == 0.9


def test_an_empty_run_writes_only_headers(tmp_path: Path):
    oligos, sets = tmp_path / "primers.tsv", tmp_path / "sets.tsv"
    assert tables.write_oligos([], oligos) == 0
    assert tables.write_sets([], sets) == 0
    assert oligos.read_text().startswith("reference\tset\trank\tprimer")
    assert sets.read_text().strip().startswith("reference\tset\trank")
