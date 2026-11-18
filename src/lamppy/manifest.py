"""Run provenance: what was run, on which inputs, with which tool versions."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from lamppy import __version__, external

_CHUNK = 1 << 20


def checksum(path: Path) -> str:
    """Return the SHA-256 digest of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def tool_versions() -> dict[str, str | None]:
    """Report the version of every external tool that is installed."""
    return {
        name: external.tool_version(name)
        for name in external.TOOLS
        if external.resolve(name) is not None
    }


def build(
    command: str,
    inputs: Iterable[Path],
    parameters: Mapping[str, Any],
) -> dict[str, Any]:
    """Assemble the provenance record of one run.

    Input checksums rather than paths are what make a result reproducible, since
    a genome file can be replaced without its name changing.
    """
    return {
        "lamppy_version": __version__,
        "command": command,
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "parameters": dict(parameters),
        "inputs": [
            {"path": str(path.resolve()), "sha256": checksum(path)} for path in inputs
        ],
        "tools": tool_versions(),
    }


def path(outdir: Path, command: str) -> Path:
    """Return where one stage records its provenance in a run directory.

    Each stage has its own file because design and validate share a directory,
    and one record overwriting the other would lose which designer made the sets.
    """
    return outdir / f"{command}.manifest.json"


def write(
    outdir: Path,
    command: str,
    inputs: Iterable[Path],
    parameters: Mapping[str, Any],
) -> Path:
    """Write the provenance record into a run directory and return its path."""
    destination = path(outdir, command)
    destination.write_text(
        json.dumps(build(command, inputs, parameters), indent=2) + "\n"
    )
    return destination
