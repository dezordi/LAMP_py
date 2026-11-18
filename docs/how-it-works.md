# How it works

```
genomes ─┬─► kmers     canonical k-mer counts, spectrum, sharing, per-genome summary
         │
         └─► targets   two-sided uniqueness gate ──► BED of candidate regions
                                                     │
                            GLAPD Single ─► par.pl ─► LAMP   (design)
                                                     │
                            thermodynamics + score + rank
                                                     │
                            blastn-short ─► configuration test   (validate)
                                                     │
                                                   report
```

## Canonical k-mers

Every k-mer is stored as the lexicographically smaller of itself and its reverse complement.
Primer binding does not care which strand a site is on, so a background match on either
strand must disqualify a region. Counts are case-insensitive — a copy is a copy — but
soft-masked windows are never *candidates*.

## The two-sided uniqueness gate

A position qualifies when its k-mer occurs at most `--max-self-copies` times (default 1) in
its own genome **and** zero times in every background genome. Runs of qualifying positions
become intervals; intervals shorter than `--min-length` are dropped.

The intra-genome half is what window BLAST against the background misses: a region repeated
five times inside its own genome and absent from every background genome passes an
inter-genome test and is useless as a target. This criterion also makes a masking chain
near-redundant, which is why masking stays optional.

Contigs are split on runs of ambiguous bases, so no interval spans an assembly gap, and each
fragment carries its offset in the source contig (`contig:start`) so every reported position
translates back to original coordinates. Inputs are never modified.

## Design

Two backends produce the same `PrimerSet` objects, so the scoring, the tables, the validation
and the report do not know which one ran.

### The primer3 backend

Candidate sites are enumerated per region under the profile's length and melting-temperature
windows, one length per start position, then assembled by the geometry of a LAMP set: F3
upstream of F2, the 5' end of F1c the required distance downstream of the 5' end of F2, the
mirror image on the B side, and a total F3-to-B3 amplicon inside its window. Only the best
arm per anchor site is kept — a 2 kb interval otherwise admits hundreds of thousands of
triples that differ by a base or two — and the assembled sets are over-generated, ranked by
melting temperature, and handed to the full scoring.

The geometric windows could not be confirmed in Eiken's documentation, so they were
calibrated against 29 sets designed by GLAPD and widened to also admit the distances quoted
for PrimerExplorer. Their provenance is in the profile's `sources`.

This backend does not check specificity. The target interval is already absent from every
background genome by construction, but no primer has been aligned anywhere, so
`lamppy validate` is the specificity check for sets designed this way.

### The GLAPD backend

All three GLAPD stages run: `Single` finds candidate primers, `par.pl` aligns them against
the background (and, with `--group`, requires them in every target-group genome) with
`bowtie`, and `LAMP` assembles the survivors into sets. Skipping `par.pl` yields sets that
were never compared with anything.

One interval is designed per run, because `par.pl` concatenates every record of its
reference and would otherwise report positions that do not exist in the genome.

Four GLAPD behaviours are guarded explicitly, each found by running it:

- **The exit status of a GLAPD stage carries no information.** `Single` returns 152 and
  `LAMP` returns 50 after completing their work normally. Every stage is therefore judged by
  the files it should have written, and a failing status is reported only when those files
  are absent.
- `par.pl` reports argument errors on standard output and still exits zero.
- `LAMP` applies `-dir` only to the files it reads, writing `-out` relative to the process
  directory, so all intervals would overwrite one file. The output path is passed absolute.
- `Single` can produce candidates for one of the two stages only, after which `LAMP` writes
  an empty file. LAMPpy checks between the stages and names the real cause.

## Scoring

Each set is measured with `primer3-py`: melting temperature, GC and end stability per
region, and hairpins and duplexes over the oligos that are actually pipetted — FIP is one
molecule, not F1c and F2 separately. Both orientations of every pair are measured, because
only a 3'-anchored duplex extends.

The penalty is a weighted sum of deviations from the [profile](parameters.md); zero means the
set satisfies every criterion. Every component is a column in `sets.tsv`, so a ranking can be
audited or recomputed without rerunning anything. Nothing is discarded by the score.

## Validation

`blastn -task blastn-short` finds every alignment of every oligo in each background, host or
database. Short queries against a whole genome produce many partial matches, so an alignment
counts as a *binding site* only when it covers at least `--min-coverage` of the oligo, has at
most `--max-mismatches` differences and reaches the oligo's 3' end — an unpaired 3' end is
not extended. Without that filter, six spurious 10 bp matches in random sequence are enough
to fake a configuration. The verdict is then about *arrangement*: does any subject sequence carry F3, F2,
F1c, B1c, B2 and B3 in that order, on the strands LAMP nomenclature implies, inside one
amplicon? Scattered partial matches are reported as `partial`; six ordered sites are
reported as `configuration`, with their positions, strands and mismatch counts as evidence.

This works on sets LAMPpy never designed — paste a published set into a two-column TSV and
screen it.
