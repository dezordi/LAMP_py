"""A single-file HTML report over the tables a run has written.

The report is assembled from the run directory alone, so it can be regenerated
after the fact and never recomputes anything. Nothing is hidden: the ranked
sets, the oligos to order, the validation evidence and the manifest all appear,
together with the statement that these are candidate assays.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path
from typing import Any

from lamppy import fasta, manifest

_SETS = "sets.tsv"
_OLIGOS = "primers.tsv"
_VALIDATION = "validation.tsv"
_STAGES = ("design", "validate")
_REPORT = "report.html"
_ORDER = ("F3", "FIP", "BIP", "B3", "LF", "LB")
_DIAGRAM_ORDER = ("F3", "F2", "F1c", "B1c", "B2", "B3")
_OLIGO_OF = {
    "F3": "F3",
    "F2": "FIP",
    "F1c": "FIP",
    "B1c": "BIP",
    "B2": "BIP",
    "B3": "B3",
}
_COLOURS = {"F3": "#6b7a8f", "FIP": "#2d6cdf", "BIP": "#d9822b", "B3": "#6b7a8f"}
_REGIONS = {
    "F3": "F3",
    "FIP": "F1c + F2",
    "BIP": "B1c + B2",
    "B3": "B3",
    "LF": "loop between F2 and F1c",
    "LB": "loop between B1c and B2",
}
_WIDTH = 760
_CAVEAT = (
    "These are candidate assays. Bioinformatic prediction is not a sufficient "
    "predictor of specificity or efficiency, and wet-lab validation is required "
    "before any diagnostic use."
)
_DESIGNERS = {
    "glapd": (
        "GLAPD, in three stages: <code>Single</code> lists candidate regions for "
        "every primer, <code>par.pl</code> drops those that align to a background "
        "genome with bowtie (and, with <code>--group</code>, those missing from a "
        "group genome), and <code>LAMP</code> combines the survivors into sets."
    ),
    "primer3": (
        "LAMPpy's built-in designer: candidate regions are enumerated within the "
        "profile's length and melting-temperature windows with primer3-py and "
        "assembled by the geometry of a LAMP set. It does not screen the "
        "background itself; specificity comes from the target interval, already "
        "absent from every background genome, and from <code>lamppy validate</code>."
    ),
}
_SET_MEANING = (
    "A LAMP set binds six regions of the target, F3, F2, F1c, B1c, B2 and B3, "
    "with four oligos: F3 and B3 bind on their own, while FIP is F1c joined to F2 "
    "and BIP is B1c joined to B2, each one molecule. The loop primers LF and LB "
    "are added when they were designed. The diagram of a set therefore shows six "
    "regions, coloured by the oligo they belong to, and its order sheet lists the "
    "four oligos to order. Every set lies inside one target interval, a region "
    "that <code>lamppy targets</code> found single-copy in this genome and absent "
    "from every background genome. A heading such as <em>Set 6 of c1:44980</em> "
    "is the set the designer numbered 6 in the interval that starts at position "
    "44980 of contig <code>c1</code>: the number identifies the set and says "
    "nothing about its quality."
)
_RANKING = (
    "Whichever designer made them, every set is measured the same way with "
    "primer3-py and compared with the design profile. Each departure from the "
    "profile is a penalty component, zero when the set meets the profile and "
    "growing with the distance from it. <code>score</code> is the weighted sum of "
    "the components and rank 1 is the lowest score. Each component is a column of "
    "the table, so the sets can be re-ranked with other weights without designing "
    "again. A low score means closer to the profile, not a working assay."
)
_STYLE = """
body { font-family: system-ui, sans-serif; margin: 2rem auto; max-width: 70rem;
       padding: 0 1rem; color: #1a1a1a; }
h1, h2 { font-weight: 600; }
table { border-collapse: collapse; width: 100%; font-size: 0.85rem;
        margin-bottom: 1.5rem; }
th, td { border-bottom: 1px solid #ddd; padding: 0.3rem 0.5rem; text-align: left; }
th { background: #f4f4f4; position: sticky; top: 0; }
code { font-family: ui-monospace, monospace; }
.caveat { background: #fff6e0; border-left: 4px solid #e0a800; padding: 0.8rem; }
.configuration { color: #a30000; font-weight: 600; }
.clear { color: #0a6b2d; }
th[title] { cursor: help; text-decoration: underline dotted; }
"""
_COLUMNS = {
    "key": "Manifest field recorded when the run was made.",
    "value": "The recorded value, as JSON.",
    "path": "Input file as it was given to the run.",
    "sha256": "SHA-256 checksum of the input; identifies the file even if it was "
    "renamed or replaced.",
    "reference": "Target interval the set was designed in, as contig:start.",
    "set": "Identifier of the primer set.",
    "rank": "Position after sorting by increasing score; 1 is the best set.",
    "contig": "Source contig the set lies on.",
    "start": "Zero-based start on the source contig.",
    "end": "Zero-based, exclusive end on the source contig.",
    "amplicon": "Length in bp from the start of F3 to the end of B3.",
    "score": "Weighted sum of the penalty_* columns; lower is better.",
    "min_conservation": "Lowest fraction of group genomes holding any one binding "
    "site of the set.",
    "worst_hairpin_dg": "Most stable hairpin among the ordered oligos, kcal/mol; "
    "more negative is worse.",
    "worst_dimer_dg": "Most stable duplex between any two oligos, kcal/mol; more "
    "negative is worse.",
    "worst_end_dimer_dg": "Most stable 3'-anchored duplex between any two oligos, "
    "kcal/mol; these can extend.",
    "primer": "Oligo name: F3, FIP, BIP, B3, and LF/LB when loop primers exist.",
    "length": "Oligo length in nucleotides.",
    "sequence": "Oligo sequence 5' to 3', as ordered.",
    "role": "Kind of source screened: background, host or database.",
    "source": "File or database the set was screened against.",
    "subject": "Sequence holding the amplifiable configuration, if any.",
    "verdict": "configuration: six sites that could amplify; partial: sites "
    "without that arrangement; clear: no site found.",
    "sites_found": "Oligos with an alignment that qualifies as a binding site.",
    "alignments": "Every local alignment BLAST returned for the set, before any "
    "threshold.",
    "configuration": "Whether the sites are arranged so they could amplify.",
    "oriented": "Whether every site of the configuration is on the strand its "
    "name implies.",
    "span": "Distance in bp from the first to the last site of the configuration.",
    "total_mismatches": "Mismatches summed over the sites of the configuration.",
    "worst_site_mismatches": "Mismatches at the least matching site of the "
    "configuration.",
    "three_prime_anchored_sites": "Binding sites whose alignment reaches the "
    "oligo's 3' end, so they can be extended.",
    "regions": "Binding regions of the oligo; FIP and BIP each join two.",
    "component": "Penalty column the weight applies to.",
    "weight": "Multiplier of the component in the score.",
    "meaning": "What the component measures.",
    "evidence": "Each site as name@position, strand +/-, and /Nmm mismatches.",
    "penalty_tm": "Degrees by which region melting temperatures fall outside "
    "their profile windows.",
    "penalty_tm_offset": "Distance of the dtm_* offsets from the profile's ideal.",
    "penalty_end_stability": "kcal/mol by which end stabilities exceed the "
    "profile limit.",
    "penalty_hairpin": "kcal/mol by which the worst hairpin is below the profile "
    "limit.",
    "penalty_dimer": "kcal/mol by which the worst duplex is below the profile "
    "limit.",
    "penalty_end_dimer": "kcal/mol by which the worst 3'-anchored duplex is below "
    "the profile limit.",
    "penalty_gc": "GC fraction outside the profile range, summed over regions.",
    "penalty_conservation": "One minus min_conservation.",
    "dtm_F2_F1c": "Melting temperature of F1c minus that of F2, in degrees.",
    "dtm_B2_B1c": "Melting temperature of B1c minus that of B2, in degrees.",
}
_VALIDATION_COLUMNS = {
    **_COLUMNS,
    "set": "Set screened, as reference/set.",
    "start": "Zero-based start of the first site on the subject.",
    "end": "Zero-based, exclusive end of the last site on the subject.",
}


def read_table(path: Path) -> list[dict[str, str]]:
    """Read a TSV table into rows, returning nothing if the file is absent."""
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def heading(field: str, descriptions: Mapping[str, str]) -> str:
    """Render a column heading, with its description as a tooltip when known."""
    description = descriptions.get(field)
    title = f' title="{escape(description)}"' if description else ""
    return f"<th{title}>{escape(field)}</th>"


def _table(
    rows: Sequence[Mapping[str, str]],
    columns: Sequence[str] | None = None,
    descriptions: Mapping[str, str] = _COLUMNS,
) -> str:
    if not rows:
        return "<p>No rows.</p>"
    fields = list(columns) if columns else list(rows[0])
    head = "".join(heading(field, descriptions) for field in fields)
    body = "".join(
        "<tr>"
        + "".join(
            f'<td class="{escape(row.get("verdict", ""))}">'
            f"{escape(str(row.get(field, '')))}</td>"
            for field in fields
        )
        + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def describe_set(row: Mapping[str, str], total: int) -> str:
    """Explain in one paragraph what a set heading refers to.

    The set number is the designer's own numbering within one target interval,
    not the rank, so both are stated to keep them from being read as the same.
    """
    reference, index = row.get("reference", ""), row.get("set", "")
    try:
        contig, offset = fasta.to_source_coordinate(reference, 0)
        interval = (
            f"the target interval that starts at position {offset} of contig "
            f"<code>{escape(contig)}</code>"
        )
    except ValueError:
        interval = f"the target interval <code>{escape(reference)}</code>"
    text = (
        f"Set {escape(index)} is the set the designer numbered {escape(index)} "
        f"in {interval}. Ranked {escape(row.get('rank', '?'))} of {total} by "
        f"score ({escape(row.get('score', '?'))}; lower is better)."
    )
    if row.get("amplicon"):
        text += (
            f" It amplifies {escape(row.get('contig', ''))}:"
            f"{escape(row.get('start', ''))}-{escape(row.get('end', ''))}, "
            f"{escape(row['amplicon'])} bp from F3 to B3."
        )
    return f"<p>{text}</p>"


def order_sheet(oligos: Sequence[Mapping[str, str]], reference: str, index: str) -> str:
    """Render the oligos of one set in the form they are ordered in."""
    rows = [
        {**row, "regions": _REGIONS[row["primer"]]}
        for row in oligos
        if row.get("reference") == reference
        and row.get("set") == index
        and row.get("primer") in _ORDER
    ]
    rows.sort(key=lambda row: _ORDER.index(row["primer"]))
    return _table(rows, ("primer", "regions", "length", "sequence"))


def diagram(oligos: Sequence[Mapping[str, str]], reference: str, index: str) -> str:
    """Draw the binding regions of one set to scale along its amplicon.

    Regions are coloured by the oligo that carries them, because FIP and BIP each
    span two regions and the order sheet lists oligos rather than regions.
    """
    sites = {
        row["primer"]: (int(row["start"]), int(row["length"]))
        for row in oligos
        if row.get("reference") == reference
        and row.get("set") == index
        and row.get("primer") in _DIAGRAM_ORDER
    }
    if len(sites) != len(_DIAGRAM_ORDER):
        return ""
    origin = min(start for start, _ in sites.values())
    span = max(start + length for start, length in sites.values()) - origin
    scale = _WIDTH / span
    boxes = "".join(
        f'<rect x="{(start - origin) * scale:.1f}" y="30" '
        f'width="{max(length * scale, 2):.1f}" height="18" '
        f'fill="{_COLOURS[_OLIGO_OF[name]]}">'
        f"<title>{escape(name)}, part of {_OLIGO_OF[name]}</title></rect>"
        f'<text x="{(start - origin) * scale:.1f}" y="24" font-size="11">'
        f"{escape(name)}</text>"
        for name, (start, length) in sites.items()
    )
    return (
        f'<svg viewBox="0 -4 {_WIDTH} 70" width="100%" role="img" '
        f'aria-label="binding sites of set {escape(index)}">'
        f'<line x1="0" y1="39" x2="{_WIDTH}" y2="39" stroke="#999" />'
        f"{boxes}"
        f'<text x="0" y="64" font-size="11">{origin}</text>'
        f'<text x="{_WIDTH}" y="64" font-size="11" text-anchor="end">'
        f"{origin + span} ({span} bp)</text></svg>"
    )


def _record(rundir: Path, command: str) -> dict[str, Any]:
    path = manifest.path(rundir, command)
    return dict(json.loads(path.read_text())) if path.exists() else {}


def overview(record: Mapping[str, Any], sets: int) -> str:
    """Say in one paragraph what was designed, from what, and with which tool."""
    if not record:
        return (
            "<p>No design manifest in this directory, so the designer is unknown.</p>"
        )
    parameters = record.get("parameters", {})
    backend = str(parameters.get("backend", ""))
    inputs = record.get("inputs", [])
    genome = Path(inputs[0]["path"]).name if inputs else "the target genome"
    commit = parameters.get("glapd_commit")
    tool = escape(backend) + (
        f" at commit <code>{escape(commit)}</code>" if commit else ""
    )
    loop = ", loop primers included" if parameters.get("loop") else ""
    intervals = parameters.get("intervals", "?")
    plural = "" if intervals == 1 else "s"
    return (
        f"<p>{sets} sets from {intervals} target interval{plural} of "
        f"<code>{escape(genome)}</code>, designed by <strong>{tool}</strong> with "
        f"the <code>{escape(str(parameters.get('profile', '?')))}</code> profile, up "
        f"to {parameters.get('num', '?')} sets per interval{loop}, by LAMPpy "
        f"{escape(str(record.get('lamppy_version', '?')))}.</p>"
        f"<p>{_DESIGNERS.get(backend, '')}</p>"
    )


def weights_table(record: Mapping[str, Any]) -> str:
    """List the weight and meaning of every penalty component the run used."""
    weights = record.get("parameters", {}).get("weights")
    if not weights:
        return ""
    return _table(
        [
            {
                "component": f"penalty_{name}",
                "weight": str(weight),
                "meaning": _COLUMNS.get(f"penalty_{name}", ""),
            }
            for name, weight in weights.items()
        ],
        ("component", "weight", "meaning"),
    )


def screening(record: Mapping[str, Any], rows: Sequence[Mapping[str, str]]) -> str:
    """Explain how the sets were screened and tally the verdicts."""
    parameters = record.get("parameters", {})
    coverage = parameters.get("min_coverage")
    rule = (
        f"covers at least {coverage:.0%} of the oligo with at most "
        f"{parameters.get('max_mismatches')} mismatches"
        if coverage is not None
        else "is close enough to the oligo"
    )
    span = parameters.get("max_span", "?")
    counts = dict.fromkeys(("configuration", "partial", "clear"), 0)
    for row in rows:
        verdict = row.get("verdict", "")
        counts[verdict] = counts.get(verdict, 0) + 1
    tally = ", ".join(f"{count} {state}" for state, count in counts.items())
    return (
        "<p>Each set was searched with <code>blastn-short</code> against every "
        f"source. An alignment counts as a binding site when it {rule} and reaches "
        "the oligo's 3' end. <code>configuration</code> means six sites in the "
        f"order and orientation of a LAMP set within {span} bp, which could "
        "amplify; <code>partial</code> means sites without that arrangement; "
        "<code>clear</code> means none were found. A clear verdict is evidence "
        f"about these sources only, not proof of specificity.</p><p>{tally}.</p>"
    )


def provenance(record: Mapping[str, Any]) -> str:
    """Render one stage's manifest: its fields, then its inputs with checksums."""
    fields = _table(
        [
            {"key": key, "value": json.dumps(value)}
            for key, value in record.items()
            if key != "inputs"
        ],
        ("key", "value"),
    )
    inputs = _table(
        [
            {"path": entry["path"], "sha256": entry["sha256"]}
            for entry in record.get("inputs", [])
        ],
        ("path", "sha256"),
    )
    return fields + inputs


def render(rundir: Path, top: int = 5) -> str:
    """Build the HTML of a report over one run directory."""
    sets = read_table(rundir / _SETS)
    oligos = read_table(rundir / _OLIGOS)
    validation = read_table(rundir / _VALIDATION)
    records = {stage: _record(rundir, stage) for stage in _STAGES}
    sections = [
        f"<h1>LAMPpy report: {escape(rundir.name)}</h1>",
        f'<p class="caveat">{_CAVEAT}</p>',
        "<h2>Summary</h2>",
        overview(records["design"], len(sets)),
        "<h2>What a set is</h2>",
        f"<p>{_SET_MEANING}</p>",
        "<h2>How the sets are ranked</h2>",
        f"<p>{_RANKING}</p>",
        weights_table(records["design"]),
        f"<h2>Ranked sets ({len(sets)})</h2>",
        "<p>Hover a column heading for its meaning. <code>rank</code>, not "
        "<code>set</code>, is the order of preference.</p>",
        _table(sets),
    ]
    for row in sets[:top]:
        reference, index = row.get("reference", ""), row.get("set", "")
        sections += [
            f"<h2>Set {escape(index)} of {escape(reference)}</h2>",
            describe_set(row, len(sets)),
            diagram(oligos, reference, index),
            order_sheet(oligos, reference, index),
        ]
    if validation:
        sections += [
            f"<h2>Validation ({len(validation)})</h2>",
            screening(records["validate"], validation),
            _table(validation, descriptions=_VALIDATION_COLUMNS),
        ]
    sections.append("<h2>Provenance</h2>")
    for stage, record in records.items():
        if record:
            sections += [f"<h3>{stage}</h3>", provenance(record)]
    body = "\n".join(sections)
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        f"<title>LAMPpy report: {escape(rundir.name)}</title>"
        f"<style>{_STYLE}</style></head><body>{body}</body></html>\n"
    )


def write(rundir: Path, destination: Path | None = None, top: int = 5) -> Path:
    """Write the HTML report of a run directory and return its path.

    Raises:
        FileNotFoundError: if the directory holds no table to report on.
    """
    if not (rundir / _SETS).exists() and not (rundir / _VALIDATION).exists():
        raise FileNotFoundError(
            f"{rundir} holds neither {_SETS} nor {_VALIDATION}; "
            "run 'lamppy design' or 'lamppy validate' first"
        )
    path = destination or rundir / _REPORT
    path.write_text(render(rundir, top))
    return path
