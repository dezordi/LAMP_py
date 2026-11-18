"""Tests for the run configuration."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lamppy import pipeline, profiles


def _config(tmp_path: Path, **overrides: object) -> Path:
    target = tmp_path / "a.fasta"
    background = tmp_path / "b.fasta"
    for path in (target, background):
        path.write_text(">c1\nACGT\n")
    document = {
        "targets": [str(target)],
        "background": [str(background)],
        "glapd": "/opt/GLAPD",
        **overrides,
    }
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump(document))
    return path


def test_glapd_is_optional_so_a_run_can_stop_at_the_targets(tmp_path: Path):
    path = _config(tmp_path)
    document = yaml.safe_load(path.read_text())
    del document["glapd"]
    path.write_text(yaml.safe_dump(document))
    assert pipeline.load(path)["glapd"] is None


def test_a_configured_glapd_path_is_read_as_a_path(tmp_path: Path):
    assert pipeline.load(_config(tmp_path))["glapd"] == Path("/opt/GLAPD")


def test_defaults_fill_in_every_stage_setting(tmp_path: Path):
    settings = pipeline.load(_config(tmp_path))
    assert settings["kmer_size"] == 31
    assert settings["profile"] == profiles.DEFAULT
    assert settings["validate"] is True
    assert settings["outdir"] == "lamppy-run"
    assert settings["group"] == ()


def test_given_values_win_over_the_defaults(tmp_path: Path):
    settings = pipeline.load(_config(tmp_path, kmer_size=21, validate=False, num=50))
    assert settings["kmer_size"] == 21
    assert settings["validate"] is False
    assert settings["num"] == 50


def test_genomes_are_read_as_paths(tmp_path: Path):
    settings = pipeline.load(_config(tmp_path))
    assert settings["targets"] == (tmp_path / "a.fasta",)
    assert settings["background"] == (tmp_path / "b.fasta",)


def test_a_single_genome_may_be_given_as_a_string(tmp_path: Path):
    path = _config(tmp_path)
    document = yaml.safe_load(path.read_text())
    document["targets"] = document["targets"][0]
    path.write_text(yaml.safe_dump(document))
    assert pipeline.load(path)["targets"] == (tmp_path / "a.fasta",)


def test_a_missing_required_key_is_named(tmp_path: Path):
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump({"targets": ["a.fasta"]}))
    with pytest.raises(ValueError, match="omits required key\\(s\\): background"):
        pipeline.load(path)


def test_an_empty_configuration_is_rejected(tmp_path: Path):
    path = tmp_path / "run.yaml"
    path.write_text("")
    with pytest.raises(ValueError, match="omits required key"):
        pipeline.load(path)


def test_a_configuration_that_is_not_a_mapping_is_rejected(tmp_path: Path):
    path = tmp_path / "run.yaml"
    path.write_text("- one\n- two\n")
    with pytest.raises(ValueError, match="not a mapping"):
        pipeline.load(path)


def test_an_absent_genome_is_reported_before_anything_runs(tmp_path: Path):
    path = _config(tmp_path, group=[str(tmp_path / "absent.fasta")])
    with pytest.raises(ValueError, match="absent.fasta"):
        pipeline.load(path)
