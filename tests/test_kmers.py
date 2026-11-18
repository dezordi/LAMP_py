"""Tests for canonical k-mer counting and sharing statistics."""

from __future__ import annotations

from pathlib import Path

import pytest

from lamppy import kmers


def test_reverse_complement_reverses_and_complements():
    assert kmers.reverse_complement("AACGT") == "ACGTT"


def test_reverse_complement_preserves_case():
    assert kmers.reverse_complement("aaCGT") == "ACGtt"


def test_canonical_picks_the_smaller_strand():
    assert kmers.canonical("TTTTT") == "AAAAA"
    assert kmers.canonical("AAAAA") == "AAAAA"


def test_canonical_collapses_both_strands_to_one_key():
    assert kmers.canonical("ACGTTGCA") == kmers.canonical("TGCAACGT")


def test_canonical_is_case_insensitive():
    assert kmers.canonical("acgtt") == kmers.canonical("ACGTT")


def test_iter_kmers_yields_every_window_with_its_offset():
    assert list(kmers.iter_kmers("AACCGG", 4)) == [
        (0, kmers.canonical("AACC")),
        (1, kmers.canonical("ACCG")),
        (2, kmers.canonical("CCGG")),
    ]


def test_iter_kmers_yields_nothing_when_the_sequence_is_shorter_than_k():
    assert list(kmers.iter_kmers("ACGT", 5)) == []


def test_count_genome_counts_a_repeated_block_twice(write_genome, random_sequence):
    block = random_sequence(100, 1)
    path = write_genome("g", {"c1": block + random_sequence(100, 2) + block})
    counts = kmers.count_genome(path, 21)
    assert counts[kmers.canonical(block[:21])] == 2


def test_count_genome_counts_both_orientations_as_one_kmer(
    write_genome, random_sequence
):
    block = random_sequence(100, 3)
    path = write_genome(
        "g", {"c1": block + random_sequence(100, 4) + kmers.reverse_complement(block)}
    )
    counts = kmers.count_genome(path, 21)
    assert counts[kmers.canonical(block[:21])] == 2


def test_count_genome_ignores_fragments_shorter_than_k(write_genome, random_sequence):
    path = write_genome("g", {"c1": random_sequence(50, 5), "c2": "ACGT"})
    assert sum(kmers.count_genome(path, 21).values()) == 30


def test_count_genomes_keys_by_file_name_not_stem(tmp_path: Path, random_sequence):
    sequence = random_sequence(60, 6)
    for name in ("g1.fa", "g1.fasta"):
        (tmp_path / name).write_text(f">c1\n{sequence}\n")
    counts = kmers.count_genomes([tmp_path / "g1.fa", tmp_path / "g1.fasta"], 21)
    assert set(counts) == {"g1.fa", "g1.fasta"}


def test_sharing_rows_omits_private_kmers(write_genome, random_sequence):
    shared = random_sequence(60, 7)
    first = write_genome("g1", {"c1": shared + random_sequence(60, 8)})
    second = write_genome("g2", {"c1": shared + random_sequence(60, 9)})
    counts = kmers.count_genomes([first, second], 21)
    rows = list(kmers.sharing_rows(counts))
    assert len(rows) == 40
    assert all(genomes == 2 for _, genomes, _, _ in rows)


def test_sharing_rows_reports_per_genome_occurrences(write_genome, random_sequence):
    shared = random_sequence(60, 10)
    first = write_genome("g1", {"c1": shared + random_sequence(60, 11) + shared})
    second = write_genome("g2", {"c1": shared + random_sequence(60, 12)})
    counts = kmers.count_genomes([first, second], 21)
    rows = {row[0]: row for row in kmers.sharing_rows(counts)}
    kmer, genomes, total, per_genome = rows[kmers.canonical(shared[:21])]
    assert (genomes, total, per_genome) == (2, 3, (2, 1))


def test_sharing_rows_honours_min_genomes(write_genome, random_sequence):
    path = write_genome("g1", {"c1": random_sequence(60, 13)})
    counts = kmers.count_genomes([path], 21)
    assert list(kmers.sharing_rows(counts, min_genomes=1)) != []
    assert list(kmers.sharing_rows(counts, min_genomes=2)) == []


def test_write_sharing_table_writes_a_header_and_one_row_per_kmer(
    tmp_path: Path, write_genome, random_sequence
):
    shared = random_sequence(60, 14)
    first = write_genome("g1", {"c1": shared})
    second = write_genome("g2", {"c1": shared})
    counts = kmers.count_genomes([first, second], 21)
    destination = tmp_path / "sharing.tsv"
    written = kmers.write_sharing_table(counts, destination)
    lines = destination.read_text().splitlines()
    assert lines[0] == "kmer\tn_genomes\ttotal_count\tg1.fasta\tg2.fasta"
    assert len(lines) == written + 1
    assert lines[1].split("\t")[1:] == ["2", "2", "1", "1"]


def test_genome_frequency_collapses_copies_within_a_genome(
    write_genome, random_sequence
):
    block = random_sequence(60, 40)
    path = write_genome("g", {"c1": block + random_sequence(60, 41) + block})
    counts = kmers.count_genomes([path], 21)
    frequency = kmers.genome_frequency(counts)
    assert counts[path.name][kmers.canonical(block[:21])] == 2
    assert frequency[kmers.canonical(block[:21])] == 1


def test_spectrum_separates_private_from_shared_kmers(write_genome, random_sequence):
    shared = random_sequence(60, 42)
    first = write_genome("g1", {"c1": shared})
    second = write_genome("g2", {"c1": shared + random_sequence(60, 43)})
    distribution = kmers.spectrum(kmers.count_genomes([first, second], 21))
    assert distribution == {1: 60, 2: 40}


def test_spectrum_of_a_clonal_set_has_no_private_kmers(write_genome, random_sequence):
    sequence = random_sequence(300, 44)
    first = write_genome("g1", {"c1": sequence})
    second = write_genome("g2", {"c1": sequence})
    distribution = kmers.spectrum(kmers.count_genomes([first, second], 21))
    assert 1 not in distribution
    assert distribution[2] == 280


def test_write_spectrum_fractions_sum_to_one(
    tmp_path: Path, write_genome, random_sequence
):
    first = write_genome("g1", {"c1": random_sequence(300, 45)})
    second = write_genome("g2", {"c1": random_sequence(300, 46)})
    counts = kmers.count_genomes([first, second], 21)
    path = tmp_path / "spectrum.tsv"
    assert kmers.write_spectrum(counts, path) == 1
    lines = path.read_text().splitlines()
    assert lines[0] == "n_genomes\tn_kmers\tfraction"
    assert sum(float(line.split("\t")[2]) for line in lines[1:]) == pytest.approx(1.0)


def test_sharing_matrix_reports_one_row_per_unordered_pair(
    write_genome, random_sequence
):
    paths = [
        write_genome(f"g{index}", {"c1": random_sequence(60, 50 + index)})
        for index in range(3)
    ]
    pairs = kmers.sharing_matrix(kmers.count_genomes(paths, 21))
    assert [(first, second) for first, second, _, _ in pairs] == [
        ("g0.fasta", "g1.fasta"),
        ("g0.fasta", "g2.fasta"),
        ("g1.fasta", "g2.fasta"),
    ]


def test_sharing_matrix_jaccard_is_one_for_identical_genomes(
    write_genome, random_sequence
):
    sequence = random_sequence(300, 47)
    first = write_genome("g1", {"c1": sequence})
    second = write_genome("g2", {"c1": sequence})
    pairs = kmers.sharing_matrix(kmers.count_genomes([first, second], 21))
    assert pairs[0][2] == 280
    assert pairs[0][3] == pytest.approx(1.0)


def test_sharing_matrix_jaccard_is_zero_for_disjoint_genomes(
    write_genome, random_sequence
):
    first = write_genome("g1", {"c1": random_sequence(300, 48)})
    second = write_genome("g2", {"c1": random_sequence(300, 49)})
    pairs = kmers.sharing_matrix(kmers.count_genomes([first, second], 21))
    assert pairs[0][2:] == (0, 0.0)


def test_summary_rows_report_copy_number_and_sharing(write_genome, random_sequence):
    block = random_sequence(60, 60)
    shared = random_sequence(60, 61)
    first = write_genome("g1", {"c1": block + shared + block})
    second = write_genome("g2", {"c1": shared})
    rows = {
        row[0]: row
        for row in kmers.summary_rows(kmers.count_genomes([first, second], 21))
    }
    _, total, distinct, private, shared_count, max_copies, multicopy = rows["g1.fasta"]
    assert total == 160
    assert max_copies == 2
    assert private + shared_count == distinct
    assert shared_count == 40
    assert multicopy == pytest.approx(80 / 160)


def test_write_summary_writes_one_row_per_genome(
    tmp_path: Path, write_genome, random_sequence
):
    first = write_genome("g1", {"c1": random_sequence(300, 62)})
    second = write_genome("g2", {"c1": random_sequence(300, 63)})
    counts = kmers.count_genomes([first, second], 21)
    path = tmp_path / "summary.tsv"
    assert kmers.write_summary(counts, path) == 2
    lines = path.read_text().splitlines()
    assert lines[0].startswith("genome\ttotal_kmers\tdistinct_kmers")
    assert lines[1].split("\t")[0] == "g1.fasta"
