# LAMPpy

LAMP primer design: find regions unique to a target genome (intra- and inter-genome), design
LAMP sets in them with GLAPD, rank them, and screen the finished sets against background
sequence. `lamppy` is a `click` CLI over `src/lamppy/`; `legacy/` holds the superseded
single-file script and is excluded from formatting, linting and tests.

## Environment

The system Python is 3.9 and the project needs 3.10+. Everything runs through the virtual
environment, which `uv` manages:

```bash
uv venv --python 3.12 .venv
uv pip install -e ".[dev,docs]" --python .venv/bin/python
```

There is no bare `python` on `PATH`. Always invoke `.venv/bin/python`, `.venv/bin/pytest`,
`.venv/bin/black`, `.venv/bin/ruff`, `.venv/bin/mypy`, `.venv/bin/mkdocs`.

## Checks

All four must pass before anything is considered done:

```bash
.venv/bin/black src tests
.venv/bin/ruff check src tests
.venv/bin/mypy src
.venv/bin/pytest -q
.venv/bin/mkdocs build --strict   # when docs or docstrings changed
```

## Conventions

- **`click`** for the CLI. `cli.py` holds decorators, argument validation and the call into a
  module — orchestration, never business logic that could live in a module.
- **One module per kind of work.** Functions of the same nature belong together:
  `kmers.py` counts, `targets.py` derives intervals, `glapd.py` drives GLAPD, `thermo.py`
  measures, `scoring.py` ranks, `tables.py` writes the final tables, `validation.py` screens.
- **`black`** is the only formatter. `ruff` is for linting only.
- **Type hints everywhere**; `mypy` runs strict.
- **No comments.** Explanations of what a function does, and why a rule exists, go in short
  objective docstrings. A docstring that states a non-obvious reason is wanted; a comment
  restating the code is not.
- **`dataclass` only for a genuine data model** (`Fragment`, `Interval`, `Primer`,
  `PrimerSet`, `Hit`, `Candidate`), never as a bag of parameters.
- **`pytest`**, one behaviour per test, the test name stating the behaviour. New behaviour
  ships with its test; GitHub Actions runs them on every pull request.
- Never reference `PROPOSAL.md` from code, docstrings or docs: that document does not ship.

## Working with external tools

Every external tool is wrapped against **observed** behaviour, not documented behaviour. Run
the tool, give it a bad argument and an empty input, see what it actually does, then guard
that. Several guards exist because a stage exits zero without doing its work:

- **Exit codes from GLAPD are meaningless**: `Single` returns 152 and `LAMP` returns 50
  after succeeding → `attempt` runs a stage, keeps the failure text, and lets the guard that
  follows decide from the files on disk.
- `par.pl` prints argument errors to stdout and exits 0 → `require_specificity_files`
  checks for the files it should have written.
- `LAMP` applies `-dir` only to the files it reads and writes `-out` relative to the process
  directory → the output path is passed absolute.
- `Single` can produce candidates for one of the two stages and exit 0, after which `LAMP`
  writes an empty file → `require_candidates` runs between them and names the real cause.

Command construction is a pure function returning `list[str]`; running it is
`external.run`. That keeps the arguments testable without the tool installed. Tools resolve
through `external.require` so a missing one gives one clear message.

Not installed on this machine: `blastn`, `makeblastdb`, `bowtie`, `dustmasker`, `trf`, `kmc`.
`docker build -t lamppy:test . && docker run --rm --entrypoint bash lamppy:test -lc '...'` is
how the real tools get exercised; Docker Desktop does not share `/private/tmp`, so generate
test data inside the container with `examples/toy_genomes.py` instead of mounting it.
Tests therefore fake `external.run` with output captured from the real tools where it was
available, and mirror the real failure modes rather than ideal ones.

## Scientific invariants

These are not style; breaking them makes the output wrong.

- **Inputs are never modified.** Everything goes to `--outdir`.
- **Coordinates always translate back** to the source contig. Fragment identifiers carry
  their offset (`contig:start`); every reported position is a source-genome position.
- **Never span an assembly gap.** Contigs are split on runs of non-ACGT bases; joining
  across one invents a junction that no primer can amplify.
- **Canonical k-mers.** A background hit on either strand disqualifies a region.
- **Uniqueness has two sides**: single-copy in its own genome *and* absent from every
  background genome. Dropping the intra-genome half admits repeats that look unique.
- **Rank, never truncate.** Nothing is discarded by the score; every component is a column.
- **No boolean "specific".** Validation reports evidence, and the docs state that wet-lab
  validation is required.
- **Soft-masking is respected**: counting is case-insensitive, but a soft-masked window is
  never a candidate target.

## Committing

Do not commit unless the user asks. The current working tree is deliberately uncommitted.
