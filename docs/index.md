# Why LAMPpy

LAMPpy designs **LAMP primer sets that are specific to one genome, or to one group of
genomes, and to nothing else you give it**. It does the part that primer design tools leave
to the user: deciding *where* on a genome an assay may be placed at all.

## The problem it solves

A LAMP set binds six sites inside about 200 bp. For an assay to be specific, that whole
window must be:

- **unique inside its own genome** — a region repeated elsewhere in the target gives
  off-target amplification from the target itself;
- **absent from every background genome** — the near neighbours the assay must not detect;
- **present in every target genome**, when the assay is meant to detect a group rather than
  one isolate;
- **conserved at the binding sites**, not merely present as a region.

LAMPpy derives those windows from canonical k-mer counts over all inputs, hands them to a
LAMP designer, then measures, scores and ranks every set it gets back, and finally screens
the finished sets for the arrangement of sites that would let a background genome amplify.

## What it does not promise

!!! warning "Candidate assays, not validated assays"
    Every output of LAMPpy is a candidate. Bioinformatic prediction is not a sufficient
    predictor of specificity or efficiency: LAMP assays fail in vitro for reasons no
    sequence analysis sees, and primer-dimer artefacts are the documented cause of LAMP
    false positives. Wet-lab validation is required before any diagnostic use.

It also does not guarantee that a target exists. Two near-clonal genomes may share every
k-mer; `lamppy kmers` tells you that in seconds, and no primer design would fix it.
