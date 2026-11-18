# Contributing

```bash
git clone https://github.com/dezordi/LAMPpy && cd LAMPpy
pip install -e ".[dev,docs]"
pytest
```

Before opening a pull request:

```bash
black .          # formatting, the only formatter
ruff check .     # linting
mypy src         # strict type checking
pytest           # tests, and they run on every pull request
sphinx-build -W -b html docs docs/_build/html  # docs, warnings as errors
```

## Conventions

- `click` for the CLI; `cli.py` holds decorators, validation and the call into a module.
- One module per kind of work; functions of the same nature stay together.
- PEP 8, type hints everywhere, `black` formatting, `ruff` for linting only.
- **No comments.** Explanations of what a function does and why a rule exists belong in
  short, objective docstrings.
- `dataclass` only where there is a genuine data model.
- Tests with `pytest`, one behaviour per test, named after the behaviour.
- A wrapper around an external tool is written against *observed* behaviour: run the tool,
  check what it actually does with a bad argument or an empty input, and guard that. Several
  guards in `glapd.py` exist because a stage exits zero without doing its work.

## Releasing to PyPI

Publishing runs from GitHub Actions (`.github/workflows/publish.yml`) through PyPI trusted
publishing, so no API token is stored anywhere.

### Once

1. Create an account on [pypi.org](https://pypi.org) and enable two-factor authentication.
2. On PyPI, open **Your account → Publishing → Add a new pending publisher** and fill in:
   project name `lamppy`, owner `dezordi`, repository `LAMPpy`, workflow `publish.yml`,
   environment `pypi`. The project is created on the first upload.
3. On GitHub, open **Settings → Environments → New environment** and name it `pypi`. Adding
   yourself as a required reviewer makes every upload wait for your approval.

### Every release

1. Set the version in `src/lamppy/__init__.py`, the only place it is defined, and move the
   `[Unreleased]` entries of `CHANGELOG.md` under `## [X.Y.Z] - YYYY-MM-DD`.
2. Run the checks above, commit and push to `master`.
3. Check the package locally:

    ```bash
    rm -rf dist && uv build && uvx twine check dist/*
    ```

4. Tag and publish a GitHub release; publishing it starts the upload:

    ```bash
    git tag vX.Y.Z && git push origin vX.Y.Z
    gh release create vX.Y.Z --title "LAMPpy X.Y.Z" --generate-notes
    ```

5. Confirm from a clean environment:

    ```bash
    pip install lamppy==X.Y.Z && lamppy --version
    ```

A version number can be uploaded only once, even after it is deleted from PyPI, so a broken
release is fixed with a new version. Versions such as `1.0.0a0` are pre-releases; while no
final release exists, `pip install lamppy` still installs the latest one.
