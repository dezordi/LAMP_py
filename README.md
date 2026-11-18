# LAMPpy

Genome-scale design of group-specific LAMP primer sets.

Given target genomes and background genomes, LAMPpy finds regions that are unique both
**intra-genome** (single-copy within their own genome) and **inter-genome** (absent from
every background genome), designs LAMP sets in them, ranks the sets by a documented score,
and screens the finished sets for background configurations that could amplify.

> **Status: not validated.** Every output is a candidate assay; wet-lab validation is
> required. See [Limitations](docs/limitations.md).

## Install

The external tools (BLAST+, bowtie, GenMap, KMC, TRF, perl) come from bioconda, not PyPI:

```bash
mamba env create -f environment.yml
mamba activate lamppy
pip install -e ".[dev]"
lamppy check-deps      # which external tools are present
```

`design` also needs a built [GLAPD](https://github.com/jiqingxiaoxi/GLAPD), which is on
neither PyPI nor bioconda. Build it once and point `GLAPD` at it (or pass `--glapd`):

```bash
git clone https://github.com/jiqingxiaoxi/GLAPD.git ~/GLAPD && make -C ~/GLAPD
export GLAPD=~/GLAPD
```

Or use the container, which has everything:

```bash
docker build -t lamppy . && docker run --rm -v "$PWD:/data" lamppy check-deps
```

## Usage

```bash
lamppy kmers     genomes/*.fasta --outdir work/ -k 31        # is there anything to target?
lamppy targets   genomes/*.fasta --outdir work/ -k 31        # unique regions, one BED per genome
lamppy prepare   genomes/*.fasta --outdir work/ [--dustmasker --trf]
lamppy design    --genome genomes/a.fasta --targets work/a.targets.bed \
                 --background genomes/b.fasta --group genomes/a2.fasta \
                 --glapd /opt/GLAPD --outdir primers/       # ranked LAMP sets
lamppy design    ... --backend primer3                      # same, with no external tool
lamppy validate  --primers primers/primers.tsv --subject genomes/b.fasta \
                 --host GRCh38.fasta --outdir primers/      # configuration test
lamppy report    --run primers/                             # one HTML file
lamppy run       --config run.yaml                          # all of the above, resumable
```

## What each stage does

**`kmers`** Read `kmers.k<k>.spectrum.tsv` first: if most k-mers are held by all *n* genomes
(`fraction` near 1 at `n_genomes = n`), the set is near-clonal with no private sequence to
target, and no primer design will fix that. The sharing matrix carries
Jaccard indices, so choosing a background set needs no separate sketching tool.

**`targets`** keeps only the regions that are single-copy in their own genome
(`--max-self-copies`) and absent from all the others. Counts are canonical, so a site on
either strand disqualifies. Output is BED6 with the copy number and the sharing genomes, a
metrics table with the length and GC of every interval, and the region sequences as FASTA —
so a run can stop here and feed any other primer designer.

**`prepare`** splits contigs on ambiguity runs so no fragment spans an assembly gap, each
identifier carrying its offset (`contig:start`). With `--dustmasker` or `--trf` it writes a
soft-masked copy — the input is never modified, and masking from both tools and from the
input itself accumulates.

**`design`** has two tool backends. `--backend primer3` designs in Python and needs no
external tool. `--backend glapd` runs all three GLAPD stages — `Single`, `par.pl`, `LAMP`
so every reported set has been aligned against the background with bowtie. Each set is 
then measured with `primer3-py` and ranked: `primers.tsv` holds the oligos with FIP/BIP
and their thermodynamics, `sets.tsv` holds one ranked row per set with every penalty
component as a column.

**`validate`** searches each oligo with `blastn-short` and tests whether any background,
host or database sequence carries all six sites in the order, orientation and spacing that
amplification requires. It works on sets LAMPpy never designed.

**`report`** states which designer and profile made the sets, explains what a set is and how
it is ranked, then gives the ranked sets, a to-scale site diagram, an order sheet and the
validation evidence in one self-contained HTML file.

Every stage writes a `<command>.manifest.json`: LAMPpy version, resolved parameters, external
tool versions and the SHA-256 of every input.

## Documentation

```bash
pip install -e ".[docs]" && mkdocs serve
```

## Development

```bash
black .          # formatting
ruff check .     # linting
mypy src         # type checking
pytest           # tests
```

## License

MIT. GLAPD (GPL-2.0) is invoked as a separate process and never vendored.
