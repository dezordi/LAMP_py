"""Tests for soft-mask detection and the unmasked-input warning."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from lamppy import fasta, masking


def test_all_uppercase_is_not_soft_masked(write_genome):
    path = write_genome("g", {"c1": "ACGTACGTACGT"})
    assert masking.detect_soft_masking(path) is False


def test_all_lowercase_is_not_soft_masked(write_genome):
    path = write_genome("g", {"c1": "acgtacgtacgt"})
    assert masking.detect_soft_masking(path) is False


def test_mixed_case_in_the_sequence_is_soft_masked(write_genome):
    path = write_genome("g", {"c1": "ACGTacgtACGT"})
    assert masking.detect_soft_masking(path) is True


def test_mixed_case_only_in_the_header_is_not_soft_masked(write_genome):
    path = write_genome("g", {"Contig_One": "acgtacgt"})
    assert masking.detect_soft_masking(path) is False


def test_lowercase_sequence_with_uppercase_ambiguity_runs_is_not_soft_masked(
    write_genome,
):
    path = write_genome("g", {"c1": "acgtacgt" + "N" * 20 + "acgtacgt"})
    assert masking.detect_soft_masking(path) is False


def test_scan_stops_at_the_byte_cap(write_genome):
    path = write_genome("g", {"c1": "ACGT" * 500 + "acgt" * 500})
    assert masking.detect_soft_masking(path, max_bytes=50) is False
    assert masking.detect_soft_masking(path) is True


def test_unmasked_genomes_trigger_exactly_one_warning_naming_both_flags(
    write_genome, caplog
):
    first = write_genome("one", {"c1": "ACGTACGT"})
    second = write_genome("two", {"c1": "ACGTACGT"})
    with caplog.at_level(logging.WARNING):
        unmasked = masking.warn_if_unmasked([first, second], False, False)
    assert unmasked == [first, second]
    assert len(caplog.records) == 1
    message = caplog.records[0].getMessage()
    assert "--dustmasker" in message
    assert "--trf" in message


def test_soft_masked_genomes_do_not_warn(write_genome, caplog):
    path = write_genome("g", {"c1": "ACGTacgtACGT"})
    with caplog.at_level(logging.WARNING):
        assert masking.warn_if_unmasked([path], False, False) == []
    assert caplog.records == []


def test_requesting_masking_suppresses_the_warning(write_genome, caplog):
    path = write_genome("g", {"c1": "ACGTACGT"})
    with caplog.at_level(logging.WARNING):
        assert masking.warn_if_unmasked([path], True, False) == []
    assert caplog.records == []


def _fake_tools(monkeypatch, behaviour):
    monkeypatch.setattr(masking.external, "require", lambda name: Path(f"/bin/{name}"))

    def run(command, cwd=None):
        behaviour(command, cwd)
        return ""

    monkeypatch.setattr(masking.external, "run", run)


def test_dustmasker_command_requests_soft_masked_fasta(tmp_path: Path):
    command = masking.dustmasker_command(
        Path("/bin/dustmasker"), tmp_path / "g.fasta", tmp_path / "out.fasta"
    )
    assert command[0] == "/bin/dustmasker"
    assert command[command.index("-outfmt") + 1] == "fasta"
    assert command[command.index("-out") + 1] == str(tmp_path / "out.fasta")


def test_trf_command_suppresses_html_and_requests_the_mask(tmp_path: Path):
    command = masking.trf_command(Path("/bin/trf"), tmp_path / "g.fasta")
    assert command[:2] == ["/bin/trf", str(tmp_path / "g.fasta")]
    assert command[2:9] == ["2", "7", "7", "80", "10", "50", "500"]
    assert "-m" in command and "-h" in command


def test_trf_output_is_named_after_the_parameters(tmp_path: Path):
    path = masking.trf_output(tmp_path / "g.fasta", tmp_path)
    assert path.name == "g.fasta.2.7.7.80.10.50.500.mask"


def test_merge_masked_lowercases_positions_reported_by_any_tool():
    merged = masking.merge_masked("ACGTACGT", "ACGTacgt", "acgtACGT")
    assert merged == "acgtacgt"


def test_merge_masked_reads_ambiguity_as_masking_and_keeps_the_base():
    assert masking.merge_masked("ACGTACGT", "ACGTNNNN") == "ACGTacgt"


def test_merge_masked_without_a_tool_returns_the_input():
    assert masking.merge_masked("ACGTacgt") == "ACGTacgt"


def test_merge_masked_rejects_a_length_mismatch():
    with pytest.raises(ValueError, match="does not match"):
        masking.merge_masked("ACGT", "ACG")


def test_soft_mask_without_a_tool_copies_the_genome(write_genome, tmp_path: Path):
    genome = write_genome("g", {"c1": "ACGTacgtACGT"})
    masked = masking.soft_mask(genome, tmp_path / "work")
    assert dict(fasta.read_records(masked)) == {"c1": "ACGTacgtACGT"}
    assert genome.read_text() != masked.read_text()


def test_soft_mask_accumulates_both_tools_and_the_input(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGTACGTacgt"})
    outdir = tmp_path / "work"

    def behaviour(command, cwd):
        outdir.mkdir(parents=True, exist_ok=True)
        if command[0].endswith("dustmasker"):
            fasta.write_records(
                [("c1", "acgtACGTACGTACGT")], Path(command[command.index("-out") + 1])
            )
        else:
            fasta.write_records(
                [("c1", "ACGTNNNNACGTACGT")], masking.trf_output(genome, cwd)
            )

    _fake_tools(monkeypatch, behaviour)
    masked = masking.soft_mask(genome, outdir, True, True)
    assert dict(fasta.read_records(masked)) == {"c1": "acgtacgtACGTacgt"}


def test_soft_mask_leaves_the_input_untouched(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT"})
    before = genome.read_bytes()

    def behaviour(command, cwd):
        (tmp_path / "work").mkdir(parents=True, exist_ok=True)
        fasta.write_records(
            [("c1", "acgtACGT")], Path(command[command.index("-out") + 1])
        )

    _fake_tools(monkeypatch, behaviour)
    masking.soft_mask(genome, tmp_path / "work", True, False)
    assert genome.read_bytes() == before


def test_run_trf_fails_when_no_mask_is_written(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT"})
    _fake_tools(monkeypatch, lambda command, cwd: None)
    with pytest.raises(RuntimeError, match="no masked sequence"):
        masking.run_trf(genome, tmp_path)


def test_run_trf_tolerates_a_nonzero_exit_when_the_mask_exists(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT"})
    monkeypatch.setattr(masking.external, "require", lambda name: Path("/bin/trf"))

    def run(command, cwd=None):
        fasta.write_records([("c1", "acgtACGT")], masking.trf_output(genome, cwd))
        raise RuntimeError("trf exited 1")

    monkeypatch.setattr(masking.external, "run", run)
    assert masking.run_trf(genome, tmp_path).exists()


def test_soft_mask_rejects_a_tool_that_renamed_a_record(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT"})

    def behaviour(command, cwd):
        fasta.write_records(
            [("other", "acgtACGT")], Path(command[command.index("-out") + 1])
        )

    _fake_tools(monkeypatch, behaviour)
    with pytest.raises(ValueError, match="no record named"):
        masking.soft_mask(genome, tmp_path / "work", True, False)


def test_run_trf_reraises_a_failure_with_no_output(
    write_genome, tmp_path: Path, monkeypatch
):
    genome = write_genome("g", {"c1": "ACGTACGT"})
    monkeypatch.setattr(masking.external, "require", lambda name: Path("/bin/trf"))

    def run(command, cwd=None):
        raise RuntimeError("trf exited 1")

    monkeypatch.setattr(masking.external, "run", run)
    with pytest.raises(RuntimeError, match="trf exited 1"):
        masking.run_trf(genome, tmp_path)
