"""Tests for the short-query BLAST wrapper."""

from __future__ import annotations

from pathlib import Path

import pytest

from lamppy import blast, fasta

TABULAR = (
    "s1|F3\tchr1\t100.000\t18\t0\t0\t1\t18\t500\t517\tplus\t1e-05\t36.2\n"
    "s1|B3\tchr1\t94.444\t18\t1\t0\t1\t18\t740\t723\tminus\t1e-03\t30.1\n"
)


def test_search_command_uses_the_short_task_and_an_explicit_format(tmp_path: Path):
    command = blast.search_command(tmp_path / "q.fasta", subject=tmp_path / "s.fasta")
    assert command[:2] == ["blastn", "-task"]
    assert command[2] == "blastn-short"
    assert command[command.index("-subject") + 1] == str(tmp_path / "s.fasta")
    assert command[command.index("-outfmt") + 1].startswith("6 qseqid sseqid")
    assert command[command.index("-dust") + 1] == "no"


def test_search_command_accepts_a_database(tmp_path: Path):
    command = blast.search_command(tmp_path / "q.fasta", database="refseq")
    assert command[command.index("-db") + 1] == "refseq"
    assert "-subject" not in command


def test_search_command_rejects_an_ambiguous_target(tmp_path: Path):
    with pytest.raises(ValueError, match="exactly one"):
        blast.search_command(tmp_path / "q.fasta")
    with pytest.raises(ValueError, match="exactly one"):
        blast.search_command(tmp_path / "q.fasta", tmp_path / "s.fasta", "refseq")


def test_hits_are_parsed_column_wise():
    first, second = blast.parse_hits(TABULAR)
    assert first.query == "s1|F3"
    assert first.subject == "chr1"
    assert first.identity == 100.0
    assert first.mismatches == 0
    assert first.bitscore == 36.2
    assert second.mismatches == 1


def test_minus_strand_coordinates_are_normalised():
    first, second = blast.parse_hits(TABULAR)
    assert (first.start, first.end, first.forward) == (499, 517, True)
    assert (second.start, second.end, second.forward) == (722, 740, False)


def test_malformed_lines_are_skipped():
    assert list(blast.parse_hits("not\ta\thit\n\n")) == []


def test_hits_group_by_subject_in_positional_order():
    extra = (
        TABULAR + "s1|F2\tchr2\t100.000\t18\t0\t0\t1\t18\t10\t27\tplus\t1e-05\t36.2\n"
    )
    grouped = blast.by_subject(blast.parse_hits(extra))
    assert sorted(grouped) == ["chr1", "chr2"]
    assert [hit.start for hit in grouped["chr1"]] == [499, 722]


def test_best_returns_the_highest_scoring_hit():
    hits = list(blast.parse_hits(TABULAR))
    assert blast.best(hits) is hits[0]
    assert blast.best([]) is None


def test_search_writes_a_query_file_and_parses_the_output(tmp_path: Path, monkeypatch):
    seen: dict[str, object] = {}

    def run(command, cwd=None):
        seen["command"] = command
        seen["queries"] = dict(
            fasta.read_records(Path(command[command.index("-query") + 1]))
        )
        return TABULAR

    monkeypatch.setattr(blast.external, "run", run)
    monkeypatch.setattr(blast.external, "require", lambda name: Path(f"/bin/{name}"))
    hits = blast.search(
        {"s1|F3": "ACGTACGTACGTACGTAC"}, tmp_path / "work", subject=tmp_path / "s.fasta"
    )
    assert len(hits) == 2
    assert seen["queries"] == {"s1|F3": "ACGTACGTACGTACGTAC"}
