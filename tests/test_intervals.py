"""Tests for the two-sided uniqueness gate, merging and BED round-tripping."""

from __future__ import annotations

from pathlib import Path

from lamppy.intervals import (
    Interval,
    drop_short,
    filter_unique,
    is_unique,
    merge,
    read_bed,
    write_bed,
)

UNIQUE = Interval("c1", 0, 300, self_copies_max=1, genomes_sharing=())
MULTI_COPY = Interval("c1", 500, 800, self_copies_max=2, genomes_sharing=())
SHARED = Interval("c1", 1000, 1300, self_copies_max=1, genomes_sharing=("other",))


def test_an_interval_unique_on_both_sides_passes():
    assert is_unique(UNIQUE) is True


def test_a_multi_copy_interval_is_rejected_despite_no_background_hit():
    assert is_unique(MULTI_COPY) is False


def test_a_shared_interval_is_rejected_despite_being_single_copy():
    assert is_unique(SHARED) is False


def test_only_the_doubly_unique_interval_reaches_design():
    assert filter_unique([UNIQUE, MULTI_COPY, SHARED]) == [UNIQUE]


def test_max_self_copies_is_a_parameter_not_a_constant():
    assert filter_unique([UNIQUE, MULTI_COPY], max_self_copies=2) == [
        UNIQUE,
        MULTI_COPY,
    ]


def test_raising_max_self_copies_never_admits_a_shared_interval():
    assert filter_unique([SHARED], max_self_copies=99) == []


def test_adjacent_intervals_merge_into_one():
    merged = merge([Interval("c1", 0, 100), Interval("c1", 100, 250)])
    assert len(merged) == 1
    assert (merged[0].start, merged[0].end) == (0, 250)


def test_intervals_separated_by_a_gap_do_not_merge():
    merged = merge([Interval("c1", 0, 100), Interval("c1", 120, 200)])
    assert len(merged) == 2


def test_a_gap_within_max_gap_merges():
    merged = merge([Interval("c1", 0, 100), Interval("c1", 120, 200)], max_gap=20)
    assert len(merged) == 1


def test_intervals_on_different_contigs_never_merge():
    merged = merge([Interval("c1", 0, 100), Interval("c2", 0, 100)])
    assert len(merged) == 2


def test_merging_takes_the_worst_metrics_so_the_gate_only_tightens():
    merged = merge(
        [
            Interval("c1", 0, 100, self_copies_max=1, genomes_sharing=("a",)),
            Interval("c1", 100, 200, self_copies_max=3, genomes_sharing=("b",)),
        ]
    )
    assert merged[0].self_copies_max == 3
    assert merged[0].genomes_sharing == ("a", "b")
    assert is_unique(merged[0]) is False


def test_merge_sorts_unordered_input():
    merged = merge([Interval("c1", 100, 200), Interval("c1", 0, 100)])
    assert len(merged) == 1
    assert (merged[0].start, merged[0].end) == (0, 200)


def test_drop_short_removes_sub_minimum_intervals():
    kept = drop_short([Interval("c1", 0, 199), Interval("c1", 300, 500)], 200)
    assert [len(interval) for interval in kept] == [200]


def test_bed_round_trips_including_metrics(tmp_path: Path):
    path = tmp_path / "targets.bed"
    intervals = [UNIQUE, MULTI_COPY, SHARED]
    assert write_bed(intervals, path) == 3
    assert read_bed(path) == intervals


def test_bed_uses_zero_based_half_open_coordinates(tmp_path: Path):
    path = tmp_path / "targets.bed"
    write_bed([Interval("c1", 10, 20)], path)
    fields = path.read_text().split("\n")[0].split("\t")
    assert fields[:3] == ["c1", "10", "20"]
    assert fields[3] == "c1:10-20"


def test_read_bed_skips_comments_and_blank_lines(tmp_path: Path):
    path = tmp_path / "targets.bed"
    path.write_text("# a comment\n\nc1\t0\t50\tc1:0-50\t1\t.\n")
    assert read_bed(path) == [Interval("c1", 0, 50)]
