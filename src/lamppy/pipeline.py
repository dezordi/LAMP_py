"""Configuration of a whole run, as read from YAML.

Only the configuration lives here: every stage of a configured run is the same
command a hand-driven run invokes, so there is one implementation of each stage
and nothing to keep in step.

A configuration without a `glapd` path is a target-discovery run: it stops after
the unique regions have been written, which is what a user wants when the primers
will be designed elsewhere. Setting `backend: primer3` instead designs without
any external tool.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import yaml

from lamppy import profiles

_REQUIRED = ("targets", "background")
_DEFAULTS: dict[str, Any] = {
    "outdir": "lamppy-run",
    "glapd": None,
    "backend": "glapd",
    "group": (),
    "host": (),
    "db": (),
    "kmer_size": 31,
    "min_length": 200,
    "max_self_copies": 1,
    "profile": profiles.DEFAULT,
    "num": 10,
    "loop": False,
    "mismatches": 2,
    "common_mismatches": 0,
    "threads": 1,
    "validate": True,
    "max_span": 500,
}
_PATH_LISTS = ("targets", "background", "group", "host")


def _as_paths(value: Any) -> tuple[Path, ...]:
    if isinstance(value, str):
        return (Path(value),)
    return tuple(Path(item) for item in value)


def load(path: Path) -> dict[str, Any]:
    """Read a run configuration, fill in the defaults and validate the inputs.

    Raises:
        ValueError: if a required key is missing or a listed genome is absent.
    """
    document = yaml.safe_load(path.read_text()) or {}
    if not isinstance(document, dict):
        raise ValueError(f"{path.name} is not a mapping of settings")
    missing = [key for key in _REQUIRED if not document.get(key)]
    if missing:
        raise ValueError(f"{path.name} omits required key(s): {', '.join(missing)}")
    settings: dict[str, Any] = {**_DEFAULTS, **document}
    for key in _PATH_LISTS:
        settings[key] = _as_paths(settings[key])
    settings["db"] = tuple(settings["db"])
    if settings["glapd"] is not None:
        settings["glapd"] = Path(settings["glapd"])
    _require_files(settings["targets"] + settings["background"] + settings["group"])
    return settings


def _require_files(paths: Sequence[Path]) -> None:
    absent = [str(path) for path in paths if not path.exists()]
    if absent:
        raise ValueError(f"configured genome(s) not found: {', '.join(absent)}")
