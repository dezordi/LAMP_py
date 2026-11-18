"""Tests for the derivation of unique target intervals."""

from __future__ import annotations

import pytest

from lamppy import kmers, targets
from lamppy.intervals import Interval

K = 21


def test_a_private_single_copy_region_is_reported(write_genome, random_sequence):
    first = write_genome("g1", {"c1": random_sequence(300, 1)})
    second = write_genome("g2", {"c1": random_sequence(300, 2)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100)
    assert [(interval.contig, interval.start, interval.end) for interval in found] == [
        ("c1", 0, 300)
    ]


def test_a_region_shared_with_the_background_is_rejected(write_genome, random_sequence):
    shared = random_sequence(300, 3)
    first = write_genome("g1", {"c1": shared + random_sequence(300, 4)})
    second = write_genome("g2", {"c1": shared + random_sequence(300, 5)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100)
    assert [(interval.start, interval.end) for interval in found] == [
        (300 - K + 1, 600)
    ]


def test_a_region_shared_on_the_reverse_strand_is_rejected(
    write_genome, random_sequence
):
    shared = random_sequence(300, 6)
    first = write_genome("g1", {"c1": shared})
    second = write_genome("g2", {"c1": kmers.reverse_complement(shared)})
    counts = kmers.count_genomes([first, second], K)
    assert targets.find_targets(first, counts, K, min_length=100) == []


def test_a_region_duplicated_inside_its_own_genome_is_rejected(
    write_genome, random_sequence
):
    duplicated = random_sequence(300, 7)
    first = write_genome(
        "g1", {"c1": duplicated, "c2": duplicated, "c3": random_sequence(300, 8)}
    )
    second = write_genome("g2", {"c1": random_sequence(300, 9)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100)
    assert [interval.contig for interval in found] == ["c3"]


def test_max_self_copies_admits_a_duplicated_region(write_genome, random_sequence):
    duplicated = random_sequence(300, 10)
    first = write_genome("g1", {"c1": duplicated, "c2": duplicated})
    second = write_genome("g2", {"c1": random_sequence(300, 11)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100, max_self_copies=2)
    assert [interval.contig for interval in found] == ["c1", "c2"]
    assert {interval.self_copies_max for interval in found} == {2}


def test_soft_masked_windows_are_skipped(write_genome, random_sequence):
    sequence = random_sequence(400, 12)
    masked = sequence[:150] + sequence[150:250].lower() + sequence[250:]
    first = write_genome("g1", {"c1": masked})
    second = write_genome("g2", {"c1": random_sequence(300, 13)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=50)
    assert [(interval.start, interval.end) for interval in found] == [
        (0, 150),
        (250, 400),
    ]


def test_intervals_use_source_coordinates_across_an_assembly_gap(
    write_genome, random_sequence
):
    sequence = random_sequence(300, 14) + "N" * 50 + random_sequence(300, 15)
    first = write_genome("g1", {"c1": sequence})
    second = write_genome("g2", {"c1": random_sequence(300, 16)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100)
    assert [(interval.start, interval.end) for interval in found] == [
        (0, 300),
        (350, 650),
    ]


def test_min_length_discards_short_intervals(write_genome, random_sequence):
    first = write_genome("g1", {"c1": random_sequence(300, 17)})
    second = write_genome("g2", {"c1": random_sequence(300, 18)})
    counts = kmers.count_genomes([first, second], K)
    assert targets.find_targets(first, counts, K, min_length=301) == []


def test_reported_intervals_carry_no_sharing_genomes(write_genome, random_sequence):
    first = write_genome("g1", {"c1": random_sequence(300, 19)})
    second = write_genome("g2", {"c1": random_sequence(300, 20)})
    counts = kmers.count_genomes([first, second], K)
    found = targets.find_targets(first, counts, K, min_length=100)
    assert all(interval.genomes_sharing == () for interval in found)


def test_find_targets_raises_when_the_genome_was_not_counted(
    write_genome, random_sequence
):
    first = write_genome("g1", {"c1": random_sequence(300, 21)})
    second = write_genome("g2", {"c1": random_sequence(300, 22)})
    counts = kmers.count_genomes([second], K)
    with pytest.raises(KeyError):
        targets.find_targets(first, counts, K)


def test_qualifying_positions_reports_the_intra_genome_copy_number():
    sequence = "ACGTTGCAACGTTGCAACGTTG"
    self_counts = {kmer: 3 for _, kmer in kmers.iter_kmers(sequence, K)}
    assert targets.qualifying_positions(sequence, K, self_counts, {}) == {}
    assert targets.qualifying_positions(
        sequence, K, self_counts, {}, max_self_copies=3
    ) == {0: 3, 1: 3}


def test_extract_cuts_interval_sequences_keeping_the_offset(write_genome):
    genome = write_genome("g1", {"c1": "AAAACCCCGGGGTTTT"})
    fragments = targets.extract(genome, [Interval("c1", 4, 8), Interval("c1", 12, 16)])
    assert [(f.identifier, f.sequence) for f in fragments] == [
        ("c1:4", "CCCC"),
        ("c1:12", "TTTT"),
    ]


def test_extract_raises_for_an_unknown_contig(write_genome):
    genome = write_genome("g1", {"c1": "ACGT"})
    with pytest.raises(KeyError):
        targets.extract(genome, [Interval("absent", 0, 4)])


def test_metrics_describe_every_interval(write_genome):
    genome = write_genome("g", {"c1": "GGCCGGCCATAT" * 10})
    found = [Interval("c1", 0, 12, 1, ()), Interval("c1", 12, 24, 2, ("b.fasta",))]
    rows = targets.metrics(genome, found)
    assert [row["length"] for row in rows] == [12, 12]
    assert rows[0]["gc"] == pytest.approx(8 / 12)
    assert rows[0]["genomes_sharing"] == "."
    assert rows[1]["genomes_sharing"] == "b.fasta"
    assert rows[1]["self_copies_max"] == 2


def test_write_metrics_writes_one_row_per_interval(write_genome, tmp_path):
    genome = write_genome("g", {"c1": "ACGT" * 100})
    found = [Interval("c1", 0, 200, 1, ())]
    path = tmp_path / "g.targets.metrics.tsv"
    assert targets.write_metrics(genome, found, path) == 1
    lines = path.read_text().splitlines()
    assert lines[0].split("\t")[:5] == ["contig", "start", "end", "length", "gc"]
    assert lines[1].split("\t")[3] == "200"
