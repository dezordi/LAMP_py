# Tutorials

## A species-specific bacterial assay

Goal: detect one species, not its congeners.

```bash
# 1. Decide the background from the data, not from intuition.
lamppy kmers refs/*.fasta --outdir work/kmers -k 31
#    sharing_matrix.tsv: Jaccard above ~0.95 means a genome belongs in the background,
#    because an assay that cannot separate it is not a species assay.

# 2. Targets: the species genome against every near neighbour.
lamppy targets refs/target.fasta refs/neighbour*.fasta --outdir work/targets -k 31

# 3. Design, requiring the assay to work for other isolates of the same species.
lamppy design --genome refs/target.fasta --targets work/targets/target.targets.bed \
              --background refs/neighbour1.fasta --background refs/neighbour2.fasta \
              --group refs/isolate2.fasta --group refs/isolate3.fasta \
              --glapd /opt/GLAPD --outdir work/design --num 50
```

`--group` does two things: `par.pl` requires every primer in every group genome, and the
binding sites are scored for conservation across the group, which enters the ranking.
Over-generate with `--num 50` and let the ranking choose; it is cheap.

```bash
# 4. Screen the finished sets against more than the design-time background.
lamppy validate --primers work/design/primers.tsv \
                --subject refs/neighbour1.fasta --db refseq_bacteria \
                --host GRCh38.fasta --outdir work/design

lamppy report --run work/design
```

Read `validation.tsv` before `sets.tsv`: a set with `verdict = configuration` against a near
neighbour is disqualified however good its score is.

## A lineage-specific viral assay

Viral genomes are small, so use a shorter k and expect fewer, shorter target intervals.

```bash
lamppy kmers lineage/*.fasta other/*.fasta --outdir work/kmers -k 18
lamppy targets lineage/patient_zero.fasta other/*.fasta --outdir work/targets \
       -k 18 --min-length 150
lamppy design --genome lineage/patient_zero.fasta \
              --targets work/targets/patient_zero.targets.bed \
              --background other/*.fasta --group lineage/*.fasta \
              --glapd /opt/GLAPD --outdir work/design --mis-c 1
```

`--mis-c 1` tolerates one mismatch when matching a primer to a group genome, which is
realistic for an RNA virus. Conservation across the lineage is the decisive column here: it
is the failure mode that broke SARS-CoV-2 LAMP assays when variants emerged.

## Validating someone else's published set

No design, no target search. Put the oligos in a TSV:

```
primer	sequence
F3	CGGTTCCCTGTACTCGAA
F2	CGAAATCTTCAAACACTACG
F1c	TGCTAATTCCTTTGTTGAGGCCG
B1c	CCTGACGGAAGCAGCATTAAG
B2	TCGAACGTAACCAAAGTCGT
B3	ACGTTGGCATTCGATTCC
```

```bash
lamppy validate --primers published.tsv --subject neighbour.fasta \
                --host GRCh38.fasta --outdir work/review
```

The configuration test needs all six regions. A set published only as FIP/BIP is screened for
sites, and the report says why the arrangement test was skipped.

## Checking that validation is working: a positive control

Screen a set against **its own target genome**. The verdict must be `configuration` with six
sites, no mismatches and `oriented = True`; anything else means the oligos, the coordinates or
the thresholds are wrong, not that the assay is specific.

```bash
lamppy validate --primers work/design/primers.tsv \
                --subject genomes/target.fasta --outdir work/control
```

```
set                  verdict        sites_found  oriented  span  total_mismatches
target_c1:44980/11   configuration  6            True      206   0
  F3@45763+/0mm F2@45793+/0mm F1c@45830-/0mm B1c@45879+/0mm B2@45905-/0mm B3@45949-/0mm
```

Do this once per project, before trusting a `clear` verdict against a background.
