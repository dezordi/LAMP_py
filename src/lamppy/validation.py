"""In-silico validation of finished LAMP sets against background sequence.

A set is not specific because its individual oligos have no good match: it is
specific because no background genome carries all six sites in an arrangement
that can amplify. The test here is therefore the arrangement, applied to sets
LAMPpy designed and to sets that arrived from a paper or a collaborator.

Not every alignment is a binding site. Short queries against a large subject
produce many partial matches that no polymerase would extend, so a hit counts
as a site only when it covers nearly the whole oligo, carries few mismatches
and reaches the oligo's 3' end.

The verdict is evidence, never a boolean: bioinformatic prediction is not a
sufficient predictor of specificity, and wet-lab validation stays mandatory.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from lamppy import blast, tables

ORDER = ("F3", "F2", "F1c", "B1c", "B2", "B3")
_COMPOSITES = ("FIP", "BIP")
_EXPECTED_FORWARD = {
    "F3": True,
    "F2": True,
    "F1c": False,
    "B1c": True,
    "B2": False,
    "B3": False,
}
_VALIDATION_FIELDS = (
    "set",
    "role",
    "source",
    "subject",
    "verdict",
    "sites_found",
    "alignments",
    "configuration",
    "oriented",
    "start",
    "end",
    "span",
    "total_mismatches",
    "worst_site_mismatches",
    "three_prime_anchored_sites",
    "evidence",
)
_MAX_MISMATCHES = 3
_MIN_COVERAGE = 0.9
CONFIGURATION = "configuration"
PARTIAL = "partial"
CLEAR = "clear"


@dataclass(frozen=True, slots=True)
class OligoSet:
    """A named set of oligos to validate, keyed by primer name."""

    label: str
    oligos: Mapping[str, str]

    @property
    def regions(self) -> dict[str, str]:
        """The six oligos the configuration test needs, when the set has them."""
        return {name: self.oligos[name] for name in ORDER if name in self.oligos}

    @property
    def testable(self) -> bool:
        """Whether the set carries all six regions in the correct nomenclature."""
        return len(self.regions) == len(ORDER)


@dataclass(frozen=True, slots=True)
class Configuration:
    """Six binding sites of one set found in amplifiable arrangement."""

    subject: str
    hits: tuple[blast.Hit, ...]

    @property
    def span(self) -> int:
        """Distance from the first to the last site."""
        return self.hits[-1].end - self.hits[0].start

    @property
    def oriented(self) -> bool:
        """Whether every site lies on the strand the nomenclature implies."""
        return all(
            hit.forward == _EXPECTED_FORWARD[hit.query.split("|")[-1]]
            for hit in self.hits
        )

    @property
    def mismatches(self) -> int:
        """Total mismatched bases over the six sites."""
        return sum(hit.mismatches for hit in self.hits)


def read_primer_table(path: Path) -> list[OligoSet]:
    """Read oligo sets from a TSV with at least primer and sequence columns.

    A set column groups the rows when present, so both a LAMPpy design table and
    a hand-typed set are accepted. Composite oligos are dropped when the regions
    they were built from are present, because a composite cannot align to a
    genome in one piece.

    Raises:
        ValueError: if the table lacks the required columns.
    """
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"primer", "sequence"} <= set(rows[0]):
        raise ValueError(f"{path.name} needs 'primer' and 'sequence' columns")
    grouped: dict[str, dict[str, str]] = {}
    for row in rows:
        label = "/".join(
            value for key in ("reference", "set") if (value := row.get(key))
        )
        grouped.setdefault(label or path.stem, {})[row["primer"]] = row["sequence"]
    return [
        OligoSet(label, _without_redundant_composites(oligos))
        for label, oligos in grouped.items()
    ]


def _without_redundant_composites(oligos: Mapping[str, str]) -> dict[str, str]:
    keep = {name: sequence for name, sequence in oligos.items()}
    if all(name in keep for name in ORDER):
        for name in _COMPOSITES:
            keep.pop(name, None)
    return keep


def plausible(
    hit: blast.Hit,
    oligos: Mapping[str, str],
    max_mismatches: int = _MAX_MISMATCHES,
    min_coverage: float = _MIN_COVERAGE,
) -> bool:
    """Whether an alignment is good enough to be a binding site.

    A site must cover nearly the whole oligo, differ from it in at most a few
    bases, and reach its 3' end, since an unpaired 3' end is not extended. The
    default thresholds are deliberately permissive about the 5' end, which
    tolerates mismatches in practice.
    """
    sequence = oligos.get(hit.query.split("|")[-1])
    if sequence is None:
        return False
    covered = hit.query_end - hit.query_start + 1
    return (
        covered >= min_coverage * len(sequence)
        and hit.mismatches + hit.gaps <= max_mismatches
        and three_prime_anchored(hit, oligos)
    )


def find_configuration(
    hits: Iterable[blast.Hit], max_span: int = 500
) -> Configuration | None:
    """Find six sites of one subject sequence arranged so they could amplify.

    The sites must occur in the order F3, F2, F1c, B1c, B2, B3 along the subject
    and within max_span bases, which is the arrangement LAMP amplification
    requires. Scattered partial matches are not a configuration.
    """
    grouped = blast.by_subject(hits)
    for subject, group in grouped.items():
        by_name: dict[str, list[blast.Hit]] = {}
        for hit in group:
            by_name.setdefault(hit.query.split("|")[-1], []).append(hit)
        chain = _chain(by_name, max_span)
        if chain is not None:
            return Configuration(subject, chain)
    return None


def _chain(
    by_name: Mapping[str, list[blast.Hit]], max_span: int
) -> tuple[blast.Hit, ...] | None:
    for first in by_name.get(ORDER[0], []):
        chosen = [first]
        for name in ORDER[1:]:
            following = next(
                (hit for hit in by_name.get(name, []) if hit.start >= chosen[-1].start),
                None,
            )
            if following is None:
                break
            chosen.append(following)
        if len(chosen) == len(ORDER) and chosen[-1].end - chosen[0].start <= max_span:
            return tuple(chosen)
    return None


def three_prime_anchored(hit: blast.Hit, oligos: Mapping[str, str]) -> bool:
    """Whether an alignment reaches the 3' end of its oligo.

    An alignment that leaves the 3' end unpaired cannot be extended, so this
    distinguishes a harmless partial match from a primeable one.
    """
    sequence = oligos.get(hit.query.split("|")[-1])
    return sequence is not None and hit.query_end >= len(sequence)


def verdict(sites: int, configuration: Configuration | None) -> str:
    """Summarise one set against one source: configuration, partial or clear."""
    if configuration is not None:
        return CONFIGURATION
    return PARTIAL if sites else CLEAR


def _query_name(label: str, primer: str) -> str:
    return f"{label}|{primer}"


def queries(primer_set: OligoSet) -> dict[str, str]:
    """Name every oligo of a set so its hits can be traced back to it."""
    return {
        _query_name(primer_set.label, name): sequence
        for name, sequence in primer_set.oligos.items()
    }


def evaluate(
    primer_set: OligoSet,
    hits: Sequence[blast.Hit],
    role: str,
    source: str,
    max_span: int = 500,
    max_mismatches: int = _MAX_MISMATCHES,
    min_coverage: float = _MIN_COVERAGE,
) -> dict[str, object]:
    """Assemble the evidence row of one set against one background source.

    Only alignments that qualify as binding sites enter the arrangement test,
    while the raw alignment count is reported alongside so a set is never
    cleared silently by a threshold.
    """
    alignments = [hit for hit in hits if hit.query.startswith(f"{primer_set.label}|")]
    own = [
        hit
        for hit in alignments
        if plausible(hit, primer_set.oligos, max_mismatches, min_coverage)
    ]
    configuration = find_configuration(own, max_span) if primer_set.testable else None
    found = {hit.query.split("|")[-1] for hit in own}
    anchored = sum(1 for hit in own if three_prime_anchored(hit, primer_set.oligos))
    row: dict[str, object] = {
        "set": primer_set.label,
        "role": role,
        "source": source,
        "verdict": verdict(len(found), configuration),
        "sites_found": len(found),
        "alignments": len(alignments),
        "configuration": configuration is not None,
        "three_prime_anchored_sites": anchored,
        "subject": "",
        "oriented": "",
        "start": "",
        "end": "",
        "span": "",
        "total_mismatches": "",
        "worst_site_mismatches": "",
        "evidence": (
            ""
            if primer_set.testable
            else f"needs all of {', '.join(ORDER)} for the configuration test"
        ),
    }
    if configuration is not None:
        row.update(
            {
                "subject": configuration.subject,
                "oriented": configuration.oriented,
                "start": configuration.hits[0].start,
                "end": configuration.hits[-1].end,
                "span": configuration.span,
                "total_mismatches": configuration.mismatches,
                "worst_site_mismatches": max(
                    hit.mismatches for hit in configuration.hits
                ),
                "evidence": " ".join(
                    f"{hit.query.split('|')[-1]}@{hit.start}"
                    f"{'+' if hit.forward else '-'}/{hit.mismatches}mm"
                    for hit in configuration.hits
                ),
            }
        )
    elif found and primer_set.testable:
        row["evidence"] = "sites without an amplifiable arrangement: " + ",".join(
            sorted(found)
        )
    return row


def write_report(rows: Sequence[dict[str, object]], path: Path) -> int:
    """Write the validation evidence table and return how many rows it holds."""
    return tables.write_rows(path, _VALIDATION_FIELDS, rows)
