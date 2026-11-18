"""The final tables of a design run: one of oligos, one of ranked sets.

Both are written as TSV with one header row, so every number behind a ranking
is auditable and recomputable without rerunning the pipeline.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from lamppy.scoring import Candidate

_OLIGO_FIELDS = (
    "reference",
    "set",
    "rank",
    "primer",
    "contig",
    "start",
    "length",
    "sequence",
    "tm",
    "gc",
    "hairpin_dg",
    "end_stability_dg",
    "conservation",
)
_SET_FIELDS = (
    "reference",
    "set",
    "rank",
    "contig",
    "start",
    "end",
    "amplicon",
    "score",
    "min_conservation",
    "worst_hairpin_dg",
    "worst_dimer_dg",
    "worst_end_dimer_dg",
)
_PENALTY_PREFIX = "penalty_"
_OFFSET_PREFIX = "dtm_"
_PRECISION = 3


def _round(value: Any) -> Any:
    return round(value, _PRECISION) if isinstance(value, float) else value


def write_rows(
    path: Path, fields: Sequence[str], rows: Sequence[dict[str, Any]]
) -> int:
    """Write rows as TSV under a fixed column order and return the row count."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(
            {key: _round(value) for key, value in row.items()} for row in rows
        )
    return len(rows)


def oligo_rows(candidates: Iterable[Candidate]) -> list[dict[str, Any]]:
    """Build one row per oligo, including the composite FIP and BIP.

    FIP and BIP carry the thermodynamics of the molecule that is pipetted, while
    the regions carry their own, because each region anneals separately.
    """
    rows: list[dict[str, Any]] = []
    for rank, candidate in enumerate(candidates, start=1):
        primer_set = candidate.primer_set
        oligos = primer_set.reaction_oligos()
        for primer in primer_set.primers:
            measured = candidate.measurements.primer(primer.name)
            rows.append(
                {
                    "reference": primer_set.reference,
                    "set": primer_set.index,
                    "rank": rank,
                    "primer": primer.name,
                    "contig": primer.contig,
                    "start": primer.start,
                    "length": len(primer),
                    "sequence": primer.sequence,
                    "tm": measured.tm if measured else "",
                    "gc": measured.gc if measured else "",
                    "hairpin_dg": measured.hairpin if measured else "",
                    "end_stability_dg": measured.end_stability if measured else "",
                    "conservation": candidate.conservation.get(primer.name, ""),
                }
            )
        for name, sequence in oligos.items():
            if primer_set.primer(name) is not None:
                continue
            anchor = primer_set.primer("F1c" if name == "FIP" else "B1c")
            rows.append(
                {
                    "reference": primer_set.reference,
                    "set": primer_set.index,
                    "rank": rank,
                    "primer": name,
                    "contig": anchor.contig if anchor else "",
                    "start": anchor.start if anchor else "",
                    "length": len(sequence),
                    "sequence": sequence,
                    "conservation": "",
                }
            )
    return rows


def set_rows(candidates: Iterable[Candidate]) -> list[dict[str, Any]]:
    """Build one row per primer set, with every penalty component as a column."""
    rows: list[dict[str, Any]] = []
    for rank, candidate in enumerate(candidates, start=1):
        primer_set = candidate.primer_set
        span = primer_set.span
        anchor = primer_set.primer("F3")
        row: dict[str, Any] = {
            "reference": primer_set.reference,
            "set": primer_set.index,
            "rank": rank,
            "contig": anchor.contig if anchor else "",
            "start": span[0] if span else "",
            "end": span[1] if span else "",
            "amplicon": span[1] - span[0] if span else "",
            "score": candidate.score.total,
            "min_conservation": candidate.weakest_site,
            "worst_hairpin_dg": candidate.measurements.worst_hairpin,
            "worst_dimer_dg": candidate.measurements.worst_dimer,
            "worst_end_dimer_dg": candidate.measurements.worst_end_dimer,
        }
        row.update(
            {
                f"{_PENALTY_PREFIX}{name}": value
                for name, value in candidate.score.components.items()
            }
        )
        row.update(
            {
                f"{_OFFSET_PREFIX}{name}": value
                for name, value in candidate.measurements.offsets.items()
            }
        )
        rows.append(row)
    return rows


def set_fields(rows: Sequence[dict[str, Any]]) -> list[str]:
    """Return the column order of the set table, penalties and offsets last."""
    extra = sorted({key for row in rows for key in row} - set(_SET_FIELDS))
    return [*_SET_FIELDS, *extra]


def write_oligos(candidates: Sequence[Candidate], path: Path) -> int:
    """Write the oligo table and return how many rows it holds."""
    return write_rows(path, _OLIGO_FIELDS, oligo_rows(candidates))


def write_sets(candidates: Sequence[Candidate], path: Path) -> int:
    """Write the ranked set table and return how many rows it holds."""
    rows = set_rows(candidates)
    return write_rows(path, set_fields(rows), rows)
