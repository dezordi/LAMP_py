# API reference

Every stage is a module that reads and writes files, so each one is usable without the CLI.
Each function below is listed with its signature and its docstring; this page is generated
from the source, so it never drifts from the code.

## Using the modules directly

The CLI is only argument handling. A whole run is a handful of calls, and the sequence below
is the one `lamppy run` performs:

```python
from pathlib import Path

from lamppy import fasta, intervals, kmers, native, profiles, scoring, tables, targets

target = Path("target.fasta")
genomes = [target, Path("near.fasta"), Path("far.fasta")]
out = Path("out")
out.mkdir(exist_ok=True)

counts = kmers.count_genomes(genomes, 31)
found = targets.find_targets(target, counts, size=31, min_length=200, max_self_copies=1)
intervals.write_bed(found, out / "target.targets.bed")
targets.write_metrics(target, found, out / "target.targets.metrics.tsv")

regions = targets.extract(target, found)
fasta.write_fragments(regions, out / "target.targets.fasta")

profile = profiles.load()
designed = [
    primer_set
    for fragment in regions
    for primer_set in native.design_fragment(fragment, profile, sets=10)
]

ranked = scoring.rank(scoring.evaluate(primer_set, profile) for primer_set in designed)
tables.write_oligos(ranked, out / "primers.tsv")
tables.write_sets(ranked, out / "sets.tsv")
```

`counts` holds every genome, the target among them: `find_targets` reads the target's own
entry as the intra-genome copy number and the remaining entries as the background, which is
how one counting pass serves both halves of the uniqueness gate.

Screening the finished sets is the two steps the `validate` command runs — one BLAST search
of every oligo, then the configuration test over its hits:

```python
from lamppy import blast, validation

sets = validation.read_primer_table(out / "primers.tsv")
oligos = {
    name: sequence
    for primer_set in sets
    for name, sequence in validation.queries(primer_set).items()
}
hits = blast.search(oligos, out / "blast", subject=Path("near.fasta"))
rows = [
    validation.evaluate(primer_set, hits, "background", "near.fasta")
    for primer_set in sets
]
validation.write_report(rows, out / "validation.tsv")
```

Three conventions hold across every module, and ignoring them is what produces wrong output:

- **Coordinates are source-genome coordinates.** A `Fragment` identifier carries its offset
  as `contig:start`, and `fasta.to_source_coordinate` translates an interval-local position
  back. Nothing reports an offset into a fragment.
- **A profile is a plain dictionary**, read once with `profiles.load` and passed down. Both
  design backends and the scoring read the same one, so a set is never measured against
  parameters other than those it was designed under.
- **Nothing is discarded by its score.** `scoring.rank` orders candidates and returns all of
  them; how many to carry forward is the caller's decision.

## Targets

```{eval-rst}
.. automodule:: lamppy.kmers
   :members:

.. automodule:: lamppy.intervals
   :members:

.. automodule:: lamppy.targets
   :members:
```

## Design

```{eval-rst}
.. automodule:: lamppy.primers
   :members:

.. automodule:: lamppy.profiles
   :members:

.. automodule:: lamppy.native
   :members:

.. automodule:: lamppy.glapd
   :members:

.. automodule:: lamppy.thermo
   :members:

.. automodule:: lamppy.scoring
   :members:

.. automodule:: lamppy.conservation
   :members:
```

## Validation and reporting

```{eval-rst}
.. automodule:: lamppy.blast
   :members:

.. automodule:: lamppy.validation
   :members:

.. automodule:: lamppy.tables
   :members:

.. automodule:: lamppy.report
   :members:
```

## Support

```{eval-rst}
.. automodule:: lamppy.fasta
   :members:

.. automodule:: lamppy.masking
   :members:

.. automodule:: lamppy.external
   :members:

.. automodule:: lamppy.manifest
   :members:

.. automodule:: lamppy.pipeline
   :members:
```
