"""Tests for the configuration test and the validation report."""

from __future__ import annotations

from pathlib import Path

import pytest

from lamppy import blast, validation

SEQUENCES = {
    "F3": "CGGTTCCCTGTACTCGAA",
    "F2": "CGAAATCTTCAAACACTACG",
    "F1c": "TGCTAATTCCTTTGTTGAGGCCG",
    "B1c": "CCTGACGGAAGCAGCATTAAG",
    "B2": "TCGAACGTAACCAAAGTCGT",
    "B3": "ACGTTGGCATTCGATTCC",
}
_STRAND = {"F3": True, "F2": True, "F1c": False, "B1c": True, "B2": False, "B3": False}


def _hit(
    primer: str,
    start: int,
    label: str = "s1",
    subject: str = "chr1",
    mismatches: int = 0,
    forward: bool | None = None,
    query_end: int | None = None,
) -> blast.Hit:
    length = len(SEQUENCES.get(primer, "A" * 20))
    oriented = _STRAND.get(primer, True) if forward is None else forward
    return blast.Hit(
        query=f"{label}|{primer}",
        subject=subject,
        identity=100.0,
        length=length,
        mismatches=mismatches,
        gaps=0,
        query_start=1,
        query_end=query_end if query_end is not None else length,
        start=start,
        end=start + length,
        strand="plus" if oriented else "minus",
        bitscore=36.0,
    )


def _configuration_hits(label: str = "s1", **kwargs: object) -> list[blast.Hit]:
    return [
        _hit(primer, 100 + offset * 50, label, **kwargs)  # type: ignore[arg-type]
        for offset, primer in enumerate(validation.ORDER)
    ]


@pytest.fixture
def primer_set() -> validation.OligoSet:
    return validation.OligoSet("s1", SEQUENCES)


def test_a_set_with_the_six_regions_is_testable(primer_set: validation.OligoSet):
    assert primer_set.testable is True
    assert list(primer_set.regions) == list(validation.ORDER)


def test_a_set_of_composites_only_is_not_testable():
    composites = validation.OligoSet("s1", {"FIP": "ACGT", "BIP": "TTGG"})
    assert composites.testable is False


def test_sites_in_order_and_span_are_a_configuration():
    found = validation.find_configuration(_configuration_hits())
    assert found is not None
    assert found.subject == "chr1"
    assert found.span == 250 + len(SEQUENCES["B3"])
    assert found.oriented is True
    assert found.mismatches == 0


def test_sites_out_of_order_are_not_a_configuration():
    hits = _configuration_hits()
    scrambled = [hits[0], hits[2], hits[1], *hits[3:]]
    scrambled[1] = _hit("F1c", 120)
    scrambled[2] = _hit("F2", 400)
    assert validation.find_configuration(scrambled) is None


def test_sites_spread_beyond_one_amplicon_are_not_a_configuration():
    hits = [
        _hit(primer, 100 + offset * 400)
        for offset, primer in enumerate(validation.ORDER)
    ]
    assert validation.find_configuration(hits) is None
    assert validation.find_configuration(hits, max_span=3000) is not None


def test_sites_scattered_over_several_contigs_are_not_a_configuration():
    hits = [
        _hit(primer, 100 + offset * 50, subject=f"chr{offset}")
        for offset, primer in enumerate(validation.ORDER)
    ]
    assert validation.find_configuration(hits) is None


def test_a_missing_site_is_not_a_configuration():
    assert validation.find_configuration(_configuration_hits()[:-1]) is None


def test_a_configuration_on_unexpected_strands_is_reported_as_such():
    found = validation.find_configuration(_configuration_hits(forward=True))
    assert found is not None
    assert found.oriented is False


def test_three_prime_anchoring_distinguishes_a_primeable_match(
    primer_set: validation.OligoSet,
):
    full = _hit("F3", 100)
    truncated = _hit("F3", 100, query_end=len(SEQUENCES["F3"]) - 3)
    assert validation.three_prime_anchored(full, primer_set.oligos) is True
    assert validation.three_prime_anchored(truncated, primer_set.oligos) is False


def test_evaluate_reports_a_configuration_with_its_evidence(
    primer_set: validation.OligoSet,
):
    row = validation.evaluate(
        primer_set, _configuration_hits(), "background", "near.fasta"
    )
    assert row["verdict"] == validation.CONFIGURATION
    assert row["sites_found"] == 6
    assert row["three_prime_anchored_sites"] == 6
    assert row["subject"] == "chr1"
    assert "F3@100+/0mm" in str(row["evidence"])


def test_evaluate_reports_partial_hits_without_an_arrangement(
    primer_set: validation.OligoSet,
):
    row = validation.evaluate(
        primer_set, [_hit("F3", 100), _hit("B2", 90)], "background", "far.fasta"
    )
    assert row["verdict"] == validation.PARTIAL
    assert row["sites_found"] == 2
    assert row["span"] == ""
    assert "without an amplifiable arrangement" in str(row["evidence"])


def test_evaluate_reports_a_clear_source(primer_set: validation.OligoSet):
    row = validation.evaluate(primer_set, [], "host", "GRCh38.fasta")
    assert row["verdict"] == validation.CLEAR
    assert row["sites_found"] == 0


def test_evaluate_ignores_hits_of_another_set(primer_set: validation.OligoSet):
    row = validation.evaluate(
        primer_set, _configuration_hits(label="s2"), "background", "near.fasta"
    )
    assert row["verdict"] == validation.CLEAR


def test_evaluate_says_why_an_untestable_set_was_not_tested():
    composites = validation.OligoSet("s1", {"FIP": "ACGT", "BIP": "TTGG"})
    row = validation.evaluate(
        composites, _configuration_hits(), "background", "b.fasta"
    )
    assert row["configuration"] is False
    assert "configuration test" in str(row["evidence"])


def test_read_primer_table_groups_by_reference_and_set(tmp_path: Path):
    path = tmp_path / "primers.tsv"
    rows = ["reference\tset\tprimer\tsequence"]
    for index in (1, 2):
        rows += [
            f"c1:0\t{index}\t{name}\t{sequence}" for name, sequence in SEQUENCES.items()
        ]
    path.write_text("\n".join(rows) + "\n")
    sets = validation.read_primer_table(path)
    assert [item.label for item in sets] == ["c1:0/1", "c1:0/2"]
    assert sets[0].oligos["F3"] == SEQUENCES["F3"]


def test_read_primer_table_accepts_a_hand_typed_set(tmp_path: Path):
    path = tmp_path / "mine.tsv"
    path.write_text(
        "primer\tsequence\n"
        + "\n".join(f"{name}\t{sequence}" for name, sequence in SEQUENCES.items())
        + "\n"
    )
    sets = validation.read_primer_table(path)
    assert len(sets) == 1
    assert sets[0].label == "mine"
    assert sets[0].testable is True


def test_read_primer_table_drops_composites_covered_by_their_regions(tmp_path: Path):
    path = tmp_path / "primers.tsv"
    rows = ["set\tprimer\tsequence"]
    rows += [f"1\t{name}\t{sequence}" for name, sequence in SEQUENCES.items()]
    rows.append("1\tFIP\t" + SEQUENCES["F1c"] + SEQUENCES["F2"])
    path.write_text("\n".join(rows) + "\n")
    assert "FIP" not in validation.read_primer_table(path)[0].oligos


def test_read_primer_table_keeps_composites_when_regions_are_absent(tmp_path: Path):
    path = tmp_path / "primers.tsv"
    path.write_text("primer\tsequence\nFIP\tACGTACGT\nBIP\tTTGGTTGG\n")
    assert set(validation.read_primer_table(path)[0].oligos) == {"FIP", "BIP"}


def test_read_primer_table_rejects_a_table_without_the_columns(tmp_path: Path):
    path = tmp_path / "wrong.tsv"
    path.write_text("name\tseq\nF3\tACGT\n")
    with pytest.raises(ValueError, match="needs 'primer' and 'sequence'"):
        validation.read_primer_table(path)


def test_write_report_holds_one_row_per_set_and_source(
    primer_set: validation.OligoSet, tmp_path: Path
):
    rows = [
        validation.evaluate(primer_set, _configuration_hits(), "background", "a.fasta"),
        validation.evaluate(primer_set, [], "host", "b.fasta"),
    ]
    path = tmp_path / "validation.tsv"
    assert validation.write_report(rows, path) == 2
    lines = path.read_text().splitlines()
    assert lines[0].split("\t")[:5] == ["set", "role", "source", "subject", "verdict"]
    assert len(lines) == 3


def test_a_full_length_near_perfect_hit_is_a_binding_site(
    primer_set: validation.OligoSet,
):
    assert validation.plausible(_hit("F3", 100), primer_set.oligos) is True
    assert (
        validation.plausible(_hit("F3", 100, mismatches=2), primer_set.oligos) is True
    )


def test_a_partial_or_mismatched_hit_is_not_a_binding_site(
    primer_set: validation.OligoSet,
):
    half = len(SEQUENCES["F3"]) // 2
    assert (
        validation.plausible(_hit("F3", 100, query_end=half), primer_set.oligos)
        is False
    )
    assert (
        validation.plausible(_hit("F3", 100, mismatches=6), primer_set.oligos) is False
    )


def test_a_hit_of_an_unknown_oligo_is_not_a_binding_site(
    primer_set: validation.OligoSet,
):
    assert validation.plausible(_hit("LF", 100), {"F3": SEQUENCES["F3"]}) is False


def test_spurious_short_alignments_do_not_make_a_configuration(
    primer_set: validation.OligoSet,
):
    spurious = [
        _hit(primer, 100 + offset * 50, query_end=8)
        for offset, primer in enumerate(validation.ORDER)
    ]
    row = validation.evaluate(primer_set, spurious, "background", "random.fasta")
    assert row["verdict"] == validation.CLEAR
    assert row["sites_found"] == 0
    assert row["alignments"] == 6


def test_the_thresholds_are_adjustable(primer_set: validation.OligoSet):
    hits = [
        _hit(primer, 100 + offset * 50, mismatches=5)
        for offset, primer in enumerate(validation.ORDER)
    ]
    strict = validation.evaluate(primer_set, hits, "background", "b.fasta")
    lenient = validation.evaluate(
        primer_set, hits, "background", "b.fasta", max_mismatches=5
    )
    assert strict["verdict"] == validation.CLEAR
    assert lenient["verdict"] == validation.CONFIGURATION
