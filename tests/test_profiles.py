"""Tests for the design-criteria profiles."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from lamppy import profiles


def test_the_default_profile_is_bundled():
    assert profiles.DEFAULT in profiles.bundled()


def test_every_bundled_profile_loads_and_documents_its_parameters():
    for name in profiles.bundled():
        profile = profiles.load(name)
        assert profile["name"] == name
        assert set(profiles.parameters(profile)) <= set(profile["sources"])


def test_the_default_profile_carries_the_eiken_windows():
    profile = profiles.load()
    assert profiles.tm_range(profile, "F1c") == (64.0, 66.0)
    assert profiles.tm_range(profile, "F3") == (59.0, 61.0)
    assert profiles.tm_range(profile, "unknown") is None


def test_unconfirmed_parameters_say_so_in_their_source():
    sources = profiles.load()["sources"]
    assert "unconfirmed" in sources["gc"]
    assert "Eiken" in sources["tm"]


def test_weights_default_to_zero_for_an_unknown_component():
    profile = profiles.load()
    assert profiles.weight(profile, "conservation") > 0.0
    assert profiles.weight(profile, "invented") == 0.0


def test_a_profile_can_be_loaded_from_a_path(tmp_path: Path):
    profile = profiles.load()
    path = tmp_path / "custom.yaml"
    path.write_text(yaml.safe_dump(profile | {"name": "custom"}))
    assert profiles.load(str(path))["name"] == "custom"


def test_an_unknown_profile_names_the_bundled_ones():
    with pytest.raises(FileNotFoundError, match=profiles.DEFAULT):
        profiles.load("no-such-profile")


def test_a_profile_missing_a_parameter_is_rejected(tmp_path: Path):
    path = tmp_path / "broken.yaml"
    path.write_text(yaml.safe_dump({"name": "broken", "parameters": {"tm": {}}}))
    with pytest.raises(ValueError, match="omits parameter"):
        profiles.load(str(path))
