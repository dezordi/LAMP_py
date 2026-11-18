# LAMPpy Documentation

**LAMPpy** is a Python CLI and package for the genome-scale design of **LAMP primer sets
that are specific to one genome, or to one group of genomes, and to nothing else you give
it**. It does the part that primer design tools leave to the user: deciding *where* on a
genome an assay may be placed at all. Its `targets` command finds those regions, `design`
builds LAMP sets in them with GLAPD or a built-in primer3-py designer, and `validate` and
`report` screen and present the finished sets.

## The problem it solves

A LAMP set binds six sites inside about 200 bp. For an assay to be specific, that whole
window must be:

- **unique inside its own genome**: a region repeated elsewhere in the target gives
  off-target amplification from the target itself;
- **absent from every background genome**: the near neighbours the assay must not detect;
- **present in every target genome**, when the assay is meant to detect a group rather than
  one isolate;
- **conserved at the binding sites**, not merely present as a region.

## Main features

- **Feasibility before design.** `lamppy kmers` counts canonical k-mers over every genome
  and reports how many are private to one genome. A near-clonal set has no private sequence,
  and that is visible in seconds, before any primer is designed. The pairwise Jaccard matrix
  also tells which genomes belong in the background.

- **Two-sided uniqueness gate.** A target interval must be single-copy in its own genome
  *and* absent from every background genome, on either strand. Dropping the intra-genome
  half would admit repeats that look unique. See [How it works](how-it-works.md).

- **Coordinates that stay on the source genome.** Contigs are split on assembly gaps and
  never joined across one, every fragment carries its offset as `contig:start`, and every
  reported position is a position on the input contig. Soft-masked sequence is never a
  target.

- **Two design backends.** `--backend glapd` runs all three GLAPD stages, including the
  `par.pl` background filter; `--backend primer3` designs in Python with primer3-py and
  needs no external tool. See [Install](install.md).

- **Ranking, not filtering.** Every set is measured the same way with primer3-py and scored
  against a versioned design profile whose parameters carry their sources. Each penalty
  component is a column, so the ranking can be re-weighted without designing again. See
  [Design parameters](parameters.md).

- **Background screening as evidence.** `lamppy validate` searches the finished sets with
  `blastn-short` against background, host and database sequence and reports whether the
  six sites ever fall into an amplifiable arrangement, with every alignment counted. It
  reports evidence, never a boolean "specific".

- **One readable report.** `lamppy report` names the designer and profile, explains what a
  set is and how it is ranked, and gives a to-scale diagram and an order sheet for each top
  set. See [Output formats](outputs.md).

- **Auditable runs.** Every stage writes a `<command>.manifest.json` with the resolved
  parameters, the SHA-256 of every input and the version of every external tool.

```{admonition} Candidate assays, not validated assays
:class: warning

Every output of LAMPpy is a candidate. Bioinformatic prediction is not a sufficient
predictor of specificity or efficiency: LAMP assays fail in vitro for reasons no sequence
analysis sees, and primer-dimer artefacts are the documented cause of LAMP false positives.
Wet-lab validation is required before any diagnostic use. See [Limitations](limitations.md).
```

LAMPpy also does not guarantee that a target exists. Two near-clonal genomes may share every
k-mer; `lamppy kmers` tells you that in seconds, and no primer design would fix it.

## Documentation contents

```{toctree}
:maxdepth: 2

install
quickstart
tutorials
how-it-works
parameters
cli
outputs
report
api
limitations
changelog
contributing
```

## Quick links

- [GitHub repository](https://github.com/dezordi/LAMPpy)
- [PyPI package](https://pypi.org/project/lamppy/)
- [Issues / bugs](https://github.com/dezordi/LAMPpy/issues)

## Cite us

If you use LAMPpy in your research, please cite it:

> Dezordi, F. Z. LAMPpy: genome-scale design of group-specific LAMP primer sets.
> https://github.com/dezordi/LAMPpy

and the method and the tools your run used:

| Tool | Used by | Reference |
|---|---|---|
| LAMP | the method itself | Notomi, T. et al. (2000). Loop-mediated isothermal amplification of DNA. *Nucleic Acids Research*, 28(12), e63. <https://doi.org/10.1093/nar/28.12.e63> |
| GLAPD | `design --backend glapd` | Jia, B. et al. (2019). GLAPD: Whole Genome Based LAMP Primer Design for a Set of Target Genomes. *Frontiers in Microbiology*, 10, 2860. <https://doi.org/10.3389/fmicb.2019.02860> |
| Bowtie | `design --backend glapd`, through GLAPD's `par.pl` | Langmead, B. et al. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. *Genome Biology*, 10(3), R25. <https://doi.org/10.1186/gb-2009-10-3-r25> |
| Primer3, through primer3-py | every `design` run, for thermodynamics and scoring; `design --backend primer3` | Untergasser, A. et al. (2012). Primer3: new capabilities and interfaces. *Nucleic Acids Research*, 40(15), e115. <https://doi.org/10.1093/nar/gks596> |
| BLAST | `validate` (`blastn`, `makeblastdb`) | Altschul, S. F. et al. (1990). Basic local alignment search tool. *Journal of Molecular Biology*, 215(3), 403–410. <https://doi.org/10.1016/S0022-2836(05)80360-2> |
| BLAST+ | `validate` and `prepare --dustmasker`, which run the BLAST+ programs | Camacho, C. et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. <https://doi.org/10.1186/1471-2105-10-421> |
| DustMasker | `prepare --dustmasker` | Morgulis, A. et al. (2006). A fast and symmetric DUST implementation to mask low-complexity DNA sequences. *Journal of Computational Biology*, 13(5), 1028–1040. <https://doi.org/10.1089/cmb.2006.13.1028> |
| Tandem Repeats Finder | `prepare --trf` | Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences. *Nucleic Acids Research*, 27(2), 573–580. <https://doi.org/10.1093/nar/27.2.573> |

LAMPpy is also built on [click](https://click.palletsprojects.com),
[PyYAML](https://pyyaml.org) and [primer3-py](https://github.com/libnano/primer3-py), which
have no paper of their own. The same list is in `CITATION.cff`, which GitHub offers as
**Cite this repository**.
