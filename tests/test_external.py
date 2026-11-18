"""Tests for external tool discovery and invocation."""

from __future__ import annotations

import sys

import pytest

from lamppy import external


def test_resolve_finds_an_installed_executable():
    assert external.resolve(sys.executable) is not None


def test_resolve_returns_none_for_a_missing_executable():
    assert external.resolve("lamppy-no-such-tool") is None


def test_require_raises_for_a_missing_executable():
    with pytest.raises(FileNotFoundError, match="check-deps"):
        external.require("lamppy-no-such-tool")


def test_run_returns_stdout():
    output = external.run([sys.executable, "-c", "print('hello')"])
    assert output.strip() == "hello"


def test_run_raises_with_stderr_on_failure():
    with pytest.raises(RuntimeError, match="boom"):
        external.run(
            [sys.executable, "-c", "import sys; sys.exit(sys.stderr.write('boom'))"]
        )


def test_run_does_not_use_a_shell():
    output = external.run([sys.executable, "-c", "print('$HOME')"])
    assert output.strip() == "$HOME"


def test_tool_version_returns_none_for_a_missing_tool():
    assert external.tool_version("lamppy-no-such-tool") is None


def test_tool_version_returns_none_for_a_tool_without_a_version_flag():
    assert external.tool_version("trf") is None


def test_genmap_probe_disables_the_update_check():
    assert external.TOOLS["genmap"] == ("--version", "--version-check", "0")
