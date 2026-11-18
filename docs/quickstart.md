# Quickstart

This runs end to end with no external tools: three toy genomes, a feasibility check and a
target search.

```bash
python examples/toy_genomes.py work/genomes
```

The script writes three 60 kb genomes: `target.fasta` and `near.fasta` share a 40 kb core,
`far.fasta` is unrelated, and `target.fasta` carries a 5 kb block twice so the intra-genome
gate has something to reject.

## 1. Is there anything to target?

```bash
lamppy kmers work/genomes/*.fasta --outdir work/kmers -k 21
```

Read `kmers.k21.spectrum.tsv` first. Its `fraction` column is the share of distinct k-mers
held by exactly `n_genomes` genomes. A large fraction at `n_genomes = 1` is private sequence,
signal to work with. If nearly everything sits at `n_genomes = 3`, the set is near-clonal:
stop here, no primer design would fix it.

`kmers.k21.summary.tsv` adds the per-genome `multicopy_fraction`, which is how much of each
genome the intra-genome gate will reject before design starts.

## 2. Where can an assay go?

```bash
lamppy targets work/genomes/*.fasta --outdir work/targets -k 21 --min-length 200
```

Every other input is the background of the genome being processed, so `target.targets.bed`
holds the regions that are single-copy in `target.fasta` **and** absent from the other two.
The duplicated block does not appear: it fails the intra-genome gate.

## 3. Design, validate, report

These need GLAPD and BLAST+ (see [Install](install.md)):

```bash
lamppy design --genome work/genomes/target.fasta \
              --targets work/targets/target.targets.bed \
              --background work/genomes/near.fasta \
              --background work/genomes/far.fasta \
              --glapd /opt/GLAPD --outdir work/design

lamppy validate --primers work/design/primers.tsv \
                --subject work/genomes/near.fasta \
                --outdir work/design

lamppy report --run work/design
```

Or all of it from one file:

```bash
lamppy run --config run.yaml
```

```yaml
targets: [work/genomes/target.fasta]
background: [work/genomes/near.fasta, work/genomes/far.fasta]
glapd: /opt/GLAPD
outdir: work/run
kmer_size: 21
validate: true
```

A stage whose output is already in place is skipped, so a run resumes where it stopped.
