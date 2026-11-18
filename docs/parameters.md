# Design parameters

Design criteria live in a YAML profile, not in code. A profile separates the numbers the
scoring uses from the citation behind each one, so no parameter can be silently invented:
the `sources` block is part of the file.

```bash
lamppy design ... --profile eiken-default     # bundled, the default
lamppy design ... --profile my-criteria.yaml  # any path
```

## The bundled default

```yaml
--8<-- "src/lamppy/profiles/eiken-default.yaml"
```

## What is confirmed and what is not

The melting-temperature windows, the F2/F1c offset and the end-stability limit are taken
from Eiken Chemical's PrimerExplorer V5 manual. The GC window, the hairpin and dimer limits
and **all weights** are LAMPpy defaults that have not been traced to a primary source; their
`sources` entries say so, and they are expected to change before 1.0.

```{admonition} Melting temperatures are not directly comparable to PrimerExplorer's
:class: note

LAMPpy computes Tm with `primer3-py` under its default salt and oligo concentrations,
which are not the conditions behind Eiken's windows. The deviation is therefore
systematic across candidates: it is fine for *ranking*, and it is not a filter — no set
is discarded for being outside a window.
```

## Reweighting

A weight of zero removes a component from the total without removing its column, so you can
rank by whatever matters for your assay:

```yaml
parameters:
  weights:
    conservation: 50.0
    gc: 0.0
```
