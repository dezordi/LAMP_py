# Changelog

All notable changes to LAMPpy are recorded here, following
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `lamppy kmers`: canonical k-mer counting with a spectrum, a pairwise sharing matrix with
  Jaccard indices, a per-genome summary and the full shared-k-mer table.
- `lamppy targets`: regions that are single-copy in their own genome and absent from every
  background genome, as BED6 with copy number and sharing per interval.
- `lamppy prepare`: splitting on ambiguity runs with coordinate-preserving identifiers, plus
  optional `--dustmasker` and `--trf` soft-masking that never modifies the input.
- `lamppy design`: all three GLAPD stages including the `par.pl` specificity filter,
  `--group` commonality, `primer3-py` thermodynamics, a reweightable penalty score and a
  ranked output with every component as a column.
- A second design backend, `--backend primer3`, which designs LAMP sets in Python with
  `primer3-py` and needs no external tool, because no installable LAMP designer exists:
  GLAPD has no package in any index, and the alternatives are web tools, a Windows-only GUI
  or abandoned. Its length, temperature and geometric windows live in the profile, calibrated
  against sets designed by GLAPD.
- `lamppy targets` writes the target regions as FASTA as well as BED, and `lamppy run` stops
  after them when no designer is configured, so the regions can feed any other tool.
- `lamppy validate`: `blastn-short` site discovery and the configuration test over
  background, host and database sequence, for LAMPpy sets and for sets from elsewhere.
- `lamppy validate` counts an alignment as a binding site only at `--min-coverage` of the
  oligo, within `--max-mismatches`, and with the 3' end aligned, which is what separates a
  real background configuration from spurious short matches.
- `lamppy report`: one self-contained HTML file that names the designer and profile, explains
  what a set is and how the ranking is computed (with the weights the run used), then gives
  the ranked sets, a to-scale site diagram coloured by oligo, an order sheet naming the
  regions each oligo joins, the validation evidence with its thresholds, and the provenance.
  Every column heading carries its meaning as a tooltip.
- `lamppy run --config`: the whole pipeline from one YAML file, resuming where it stopped.
- `lamppy check-deps`: external tools and their versions.
- `<command>.manifest.json` per stage: input checksums, resolved parameters, tool versions,
  platform and timestamp. One file per stage, because `design` and `validate` share a
  directory. The design manifest records the backend, the GLAPD commit and the score weights.
- Versioned design-criteria profiles in YAML, each parameter carrying its source.
- A per-interval metrics table from `lamppy targets` with length, GC and copy number.
- A Dockerfile pinning every external tool and building GLAPD from a pinned commit.

### Changed

- The pipeline is a `click` CLI and an importable package, replacing the single script.
- Uniqueness is decided per interval over canonical k-mer counts, intra-genome and
  inter-genome, in place of window BLAST against the background only.
- Primer sets are ranked rather than truncated to the first the designer found.

### Fixed

- GLAPD stages are judged by the files they wrote rather than by their exit status:
  `Single` returns 152 and `LAMP` returns 50 after succeeding, which made every design run
  fail against the real tool.
- `par.pl` is no longer skipped, so reported sets have been compared against the background.
- GLAPD stages that exit zero without doing their work are detected and named: `par.pl`
  argument errors on standard output, `LAMP` ignoring `-dir` for `-out`, and `Single`
  producing candidates for one stage only.
- Inputs are never modified and never concatenated across assembly gaps.
