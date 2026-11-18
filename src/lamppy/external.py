"""Discovery and invocation of the external tools the pipeline delegates to."""

from __future__ import annotations

import logging
import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

logger = logging.getLogger(__name__)

TOOLS: dict[str, tuple[str, ...] | None] = {
    "blastn": ("-version",),
    "makeblastdb": ("-version",),
    "bowtie": ("--version",),
    "perl": ("--version",),
    "bowtie-build": ("--version",),
    "dustmasker": ("-version",),
    "trf": None,
}


def resolve(executable: str) -> Path | None:
    """Return the path of an executable on PATH, or None if absent."""
    found = shutil.which(executable)
    return Path(found) if found else None


def require(executable: str) -> Path:
    """Return the path of an executable on PATH.

    Raises:
        FileNotFoundError: if the executable is not installed.
    """
    found = resolve(executable)
    if found is None:
        raise FileNotFoundError(
            f"{executable} is not on PATH; run 'lamppy check-deps' for the full list"
        )
    return found


def run(command: Sequence[str], cwd: Path | None = None) -> str:
    """Run a command without a shell and return its stdout.

    Raises:
        RuntimeError: if the command exits non-zero.
    """
    logger.debug("running %s", " ".join(command))
    completed = subprocess.run(
        [str(part) for part in command],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{command[0]} failed with code {completed.returncode}: "
            f"{completed.stderr.strip()}"
        )
    return completed.stdout


def tool_version(executable: str) -> str | None:
    """Return the first line of a tool's version output, or None if unavailable."""
    arguments = TOOLS.get(executable)
    if arguments is None or resolve(executable) is None:
        return None
    try:
        output = run([executable, *arguments])
    except RuntimeError:
        return None
    return next((line.strip() for line in output.splitlines() if line.strip()), None)
