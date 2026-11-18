# Limitations and wet-lab validation

!!! danger "Required reading"
    LAMPpy outputs **candidate** assays. Bioinformatic prediction is not a sufficient
    predictor of specificity or efficiency. Wet-lab validation — inclusivity panel,
    exclusivity panel, limit of detection, and no-template controls — is required before any
    diagnostic use.

## Known limits of the method

- **Background is only what you give it.** A set is specific against the genomes and
  databases it was screened against, nothing more. The design-time background should be the
  near neighbours; `validate` is where a large database belongs.
- **Closely related targets are enriched for accessory and repetitive content.** The regions
  unique to one of two near-identical genomes tend to be mobile elements and plasmids. The
  intra-genome gate handles the repetitive half; the accessory half is a biology problem no
  masking step can fix. Screen final sets against a large database.
- **Exact-match conservation.** `--group` scores a binding site as conserved only when the
  whole oligo matches exactly. That is strict by design, and it is not an alignment: a site
  with one tolerated mismatch at the 5' end counts as absent.
- **Melting temperatures use primer3 defaults**, not LAMP reaction conditions. Fine for
  ranking, not an absolute prediction. See [Design parameters](parameters.md).
- **k-mer counts are held in memory**, which is adequate for viral and bacterial genome sets
  (around 0.6 GB per 5 Mbp genome at k=31). Eukaryotic inputs would need KMC3 or unikmer.
- **No workflow engine.** Stages skip when their output exists; they do not verify that the
  inputs behind an existing output are unchanged.
