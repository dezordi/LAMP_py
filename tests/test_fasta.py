"""Tests for FASTA reading, fragment splitting and coordinate translation."""

from __future__ import annotations

from pathlib import Path

import pytest

from lamppy import fasta


def test_read_records_joins_wrapped_lines_and_truncates_identifier(write_genome):
    path = write_genome("g", {"contig1": "ACGT" * 50})
    records = list(fasta.read_records(path))
    assert records == [("contig1", "ACGT" * 50)]


def test_read_records_preserves_case(write_genome):
    path = write_genome("g", {"contig1": "ACGTacgtACGT"})
    assert list(fasta.read_records(path))[0][1] == "ACGTacgtACGT"


def test_read_records_handles_multiple_records(write_genome):
    path = write_genome("g", {"a": "ACGTACGT", "b": "TTTTGGGG"})
    assert list(fasta.read_records(path)) == [("a", "ACGTACGT"), ("b", "TTTTGGGG")]


def test_split_covers_the_last_base_of_the_contig():
    sequence = "ACGT" * 100
    fragments = fasta.split_at_ambiguous("c", sequence)
    assert len(fragments) == 1
    assert fragments[0].start == 0
    assert len(fragments[0]) == len(sequence)
    assert fragments[0].sequence.endswith(sequence[-10:])


def test_split_never_spans_an_ambiguity_run():
    sequence = "AAAACCCC" + "N" * 10 + "GGGGTTTT"
    fragments = fasta.split_at_ambiguous("c", sequence)
    assert [(f.start, f.sequence) for f in fragments] == [
        (0, "AAAACCCC"),
        (18, "GGGGTTTT"),
    ]


def test_split_breaks_on_iupac_ambiguity_codes():
    fragments = fasta.split_at_ambiguous("c", "AAAARCCCC")
    assert [f.sequence for f in fragments] == ["AAAA", "CCCC"]


def test_split_keeps_soft_masked_bases():
    fragments = fasta.split_at_ambiguous("c", "AAAAacgtTTTT")
    assert [f.sequence for f in fragments] == ["AAAAacgtTTTT"]


def test_split_discards_fragments_below_min_length():
    sequence = "AAA" + "N" + "CCCCCCCC"
    fragments = fasta.split_at_ambiguous("c", sequence, min_length=5)
    assert [f.sequence for f in fragments] == ["CCCCCCCC"]


def test_split_of_an_all_ambiguous_sequence_is_empty():
    assert fasta.split_at_ambiguous("c", "NNNNNN") == []


def test_fragment_position_round_trips_to_the_source_contig():
    sequence = "AAAACCCC" + "N" * 10 + "GGGGTTTT"
    fragment = fasta.split_at_ambiguous("contig1", sequence)[1]
    contig, position = fasta.to_source_coordinate(fragment.identifier, 2)
    assert (contig, position) == ("contig1", 20)
    assert sequence[position] == fragment.sequence[2]


def test_to_source_coordinate_rejects_a_plain_identifier():
    with pytest.raises(ValueError, match="not a fragment identifier"):
        fasta.to_source_coordinate("contig1", 0)


def test_to_source_coordinate_tolerates_colons_in_the_contig_name():
    assert fasta.to_source_coordinate("ref:v1:300", 5) == ("ref:v1", 305)


def test_write_fragments_round_trips(tmp_path: Path, write_genome):
    path = write_genome("g", {"c1": "ACGT" * 30 + "N" * 5 + "TTTT" * 30})
    destination = tmp_path / "out.fasta"
    written = fasta.write_fragments(fasta.iter_fragments(path), destination)
    assert written == 2
    records = list(fasta.read_records(destination))
    assert [name for name, _ in records] == ["c1:0", "c1:125"]
    assert records[0][1] == "ACGT" * 30


def test_inputs_are_byte_identical_after_processing(write_genome):
    path = write_genome("g", {"c1": "ACGTN" * 40})
    before = path.read_bytes()
    list(fasta.iter_fragments(path))
    assert path.read_bytes() == before
