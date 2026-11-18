# Output formats

Every stage writes plain text into `--outdir` and never touches its inputs.

## `kmers`

| File | Columns |
|---|---|
| `kmers.k<k>.spectrum.tsv` | `n_genomes`, `n_kmers`, `fraction` |
| `kmers.k<k>.sharing_matrix.tsv` | `genome_a`, `genome_b`, `shared_kmers`, `jaccard` |
| `kmers.k<k>.summary.tsv` | `genome`, `total_kmers`, `distinct_kmers`, `private_kmers`, `shared_kmers`, `max_self_copies`, `multicopy_fraction` |
| `kmers.k<k>.sharing.tsv` | `kmer`, `n_genomes`, `total_count`, then one occurrence-count column per genome |

The first three are kilobytes and always written. The fourth is one row per shared k-mer and
can be large; `--min-genomes` bounds it.

The Jaccard index in the sharing matrix is the neighbour distance that decides which genomes
belong in a background set, so no separate sketching tool is needed.

## `targets`

BED6 per genome, `<genome>.targets.bed`: `contig`, `start`, `end`, `name`, then column 5 is
the intra-genome copy number of the interval and column 6 lists the genomes that share it
(`.` when none do). Coordinates are source-contig coordinates, zero-based half-open.

`<genome>.targets.metrics.tsv` carries one row per interval for review before any design
runs: `contig`, `start`, `end`, `length`, `gc`, `self_copies_max`, `genomes_sharing`. GC is
the metric that decides feasibility — a region too AT-rich for the melting temperatures a
LAMP set needs yields no sets however unique it is.

## `design`

`primers.tsv` — one row per oligo:

`reference`, `set`, `rank`, `primer`, `contig`, `start`, `length`, `sequence`, `tm`, `gc`,
`hairpin_dg`, `end_stability_dg`, `conservation`

The six GLAPD regions appear with their own thermodynamics; **FIP** (F1c+F2) and **BIP**
(B1c+B2) appear as the oligos you order, with the start of their first half.

`sets.tsv` — one row per set, ranked, with every number behind the ranking:

`reference`, `set`, `rank`, `contig`, `start`, `end`, `amplicon`, `score`,
`min_conservation`, `worst_hairpin_dg`, `worst_dimer_dg`, `worst_end_dimer_dg`,
`dtm_F2_F1c`, `dtm_B2_B1c`, and one `penalty_*` column per score component.

Free energies are kcal/mol. `rank` 1 is the lowest penalty.

## `validate`

`validation.tsv` — one row per set and source:

`set`, `role` (`background`, `host`, `database`), `source`, `subject`, `verdict`,
`sites_found`, `alignments`, `configuration`, `oriented`, `start`, `end`, `span`,
`total_mismatches`, `worst_site_mismatches`, `three_prime_anchored_sites`, `evidence`

`alignments` is every local alignment BLAST returned; `sites_found` counts only those that
qualify as binding sites (`--min-coverage` of the oligo, at most `--max-mismatches`
differences, 3' end aligned). Both are reported so no set is cleared silently by a
threshold.

`verdict` is `configuration` (six sites arranged so they could amplify — a real risk),
`partial` (sites present, no amplifiable arrangement) or `clear` (nothing found). `evidence`
spells out each site as `name@position±/Nmm`.

## `report`

`report.html`, one self-contained file:

- a summary naming the designer (with the GLAPD commit when known), the profile and the genome;
- what a set is: six binding regions carried by four oligos, FIP joining F1c to F2 and BIP
  joining B1c to B2, inside one target interval, numbered by the designer;
- how the ranking works, with the weight of every penalty component the run used;
- the ranked sets, then for the top sets a to-scale diagram with regions coloured by the oligo
  that carries them and an order sheet naming the regions of each oligo;
- the validation evidence, with the thresholds it was judged by and a tally of verdicts;
- the provenance of each stage.

Hovering a column heading shows what the column means.

## `<command>.manifest.json`

Written by every stage, one file each (`design.manifest.json`, `validate.manifest.json`,
...), so that `design` and `validate` can share a directory without one overwriting the
other: LAMPpy version, command, UTC timestamp, Python and platform, the
resolved parameters, the SHA-256 of every input and the version string of every external
tool found. Checksums rather than paths, because a FASTA can be replaced without its name
changing.
