# Legacy script

`script_v_0_4_5.py` is the original single-file pipeline, kept for reference and for
reproducing earlier results. **It is not maintained and should not be used for new work.**

Use `lamppy` instead. The differences that matter:

| Legacy script | LAMPpy |
|---|---|
| k-mer windows compared by BLAST against the background only | canonical k-mer counts, unique intra-genome **and** inter-genome |
| GLAPD `Single` and `LAMP`, no `par.pl` | all three GLAPD stages, so sets are compared against the background |
| first `-num` sets reported | every set measured, scored and ranked |
| no validation of finished sets | `lamppy validate` configuration test |
| inputs rewritten in place, `os.chdir` into the GLAPD directory | inputs read-only, everything written to `--outdir` |
| no record of what was run | `manifest.json` with checksums and tool versions |

Equivalent of a legacy run:

```bash
lamppy targets genomes/*.fasta --outdir work/ -k 21
lamppy design --genome genomes/a.fasta --targets work/a.targets.bed \
              --background genomes/b.fasta --glapd /opt/GLAPD --outdir primers/
```
