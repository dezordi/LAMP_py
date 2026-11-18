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
mkdocs build -s  # docs, strict
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
