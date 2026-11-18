"""Versioned design-criteria profiles, loaded from YAML.

A profile separates the numbers the scoring uses from the citation behind each
one, so a parameter can never be silently invented: the sources block is part
of the file and is rendered in the documentation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

_DIRECTORY = Path(__file__).parent / "profiles"
_SUFFIX = ".yaml"
DEFAULT = "eiken-default"
_REQUIRED = ("tm", "tm_offset", "end_stability_dg", "dimer_dg", "gc", "weights")


def bundled() -> list[str]:
    """List the profiles shipped with LAMPpy."""
    return sorted(path.stem for path in _DIRECTORY.glob(f"*{_SUFFIX}"))


def load(name: str = DEFAULT) -> dict[str, Any]:
    """Load a bundled profile by name, or any profile by path.

    Raises:
        FileNotFoundError: if neither a bundled profile nor a file matches.
        ValueError: if the profile omits a required parameter.
    """
    path = _DIRECTORY / f"{name}{_SUFFIX}"
    if not path.exists():
        path = Path(name)
    if not path.exists():
        raise FileNotFoundError(
            f"unknown profile {name!r}; bundled profiles: {', '.join(bundled())}"
        )
    profile: dict[str, Any] = yaml.safe_load(path.read_text())
    parameters = profile.get("parameters", {})
    missing = [key for key in _REQUIRED if key not in parameters]
    if missing:
        raise ValueError(f"{path.name} omits parameter(s): {', '.join(missing)}")
    return profile


def parameters(profile: dict[str, Any]) -> dict[str, Any]:
    """Return the numeric parameters of a profile."""
    values: dict[str, Any] = profile["parameters"]
    return values


def weight(profile: dict[str, Any], component: str) -> float:
    """Return the weight of one score component, or zero if it has none."""
    return float(parameters(profile)["weights"].get(component, 0.0))


def tm_range(profile: dict[str, Any], primer: str) -> tuple[float, float] | None:
    """Return the melting-temperature window of one primer region."""
    window = parameters(profile)["tm"].get(primer)
    return (float(window[0]), float(window[1])) if window else None
