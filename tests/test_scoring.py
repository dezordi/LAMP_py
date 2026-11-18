"""Tests for the penalty score and the ranking it drives."""

from __future__ import annotations

from pathlib import Path

from lamppy import glapd, profiles, scoring, thermo

PRIMERS = {
    "F3": "CGGTTCCCTGTACTCGAA",
    "F2": "CGAAATCTTCAAACACTACG",
    "F1c": "TGCTAATTCCTTTGTTGAGGCCG",
    "B1c": "CCTGACGGAAGCAGCATTAAG",
    "B2": "TCGAACGTAACCAAAGTCGT",
    "B3": "ACGTTGGCATTCGATTCC",
}


def _primer_set(
    index: int = 1, sequences: dict[str, str] | None = None
) -> glapd.PrimerSet:
    oligos = sequences or PRIMERS
    return glapd.PrimerSet(
        reference="c1:0",
        index=index,
        primers=tuple(
            glapd.Primer(name, "c1", 100 + offset * 40, sequence)
            for offset, (name, sequence) in enumerate(oligos.items())
        ),
    )


def test_every_component_is_reported():
    profile = profiles.load()
    parts = scoring.components(thermo.evaluate_set(PRIMERS), profile)
    assert set(parts) == {
        "tm",
        "tm_offset",
        "end_stability",
        "hairpin",
        "dimer",
        "end_dimer",
        "gc",
        "conservation",
    }
    assert all(value >= 0.0 for value in parts.values())


def test_a_component_is_zero_when_the_profile_is_satisfied(tmp_path: Path):
    profile = profiles.load()
    profile["parameters"]["gc"] = [0.0, 1.0]
    parts = scoring.components(thermo.evaluate_set(PRIMERS), profile)
    assert parts["gc"] == 0.0


def test_a_tm_outside_the_window_is_penalised_by_its_distance():
    profile = profiles.load()
    profile["parameters"]["tm"] = {"F3": [0.0, 10.0]}
    measured = thermo.evaluate_set({"F3": PRIMERS["F3"]})
    parts = scoring.components(measured, profile)
    assert parts["tm"] == measured.primers[0].tm - 10.0


def test_conservation_enters_the_penalty_as_a_shortfall():
    profile = profiles.load()
    measured = thermo.evaluate_set(PRIMERS)
    assert scoring.components(measured, profile, 1.0)["conservation"] == 0.0
    assert scoring.components(measured, profile, 0.25)["conservation"] == 0.75


def test_the_total_is_the_weighted_sum_of_the_components():
    profile = profiles.load()
    measured = thermo.evaluate_set(PRIMERS)
    result = scoring.score(measured, profile, 0.5)
    expected = sum(
        profiles.weight(profile, name) * value
        for name, value in result.components.items()
    )
    assert result.total == expected


def test_a_weight_of_zero_removes_a_component_from_the_total():
    profile = profiles.load()
    measured = thermo.evaluate_set(PRIMERS)
    with_conservation = scoring.score(measured, profile, 0.0).total
    profile["parameters"]["weights"]["conservation"] = 0.0
    without = scoring.score(measured, profile, 0.0).total
    assert without < with_conservation


def test_evaluate_measures_duplexes_over_the_composite_oligos():
    candidate = scoring.evaluate(_primer_set(), profiles.load())
    regions_only = thermo.evaluate_set(PRIMERS)
    assert candidate.measurements.worst_dimer != regions_only.worst_dimer
    assert candidate.weakest_site == 1.0


def test_evaluate_takes_the_weakest_site_as_the_set_conservation():
    candidate = scoring.evaluate(_primer_set(), profiles.load(), {"F3": 1.0, "B3": 0.4})
    assert candidate.weakest_site == 0.4
    assert candidate.score.components["conservation"] == 0.6


def test_rank_puts_the_lowest_penalty_first():
    profile = profiles.load()
    poor = dict(PRIMERS) | {"F3": "AAAAAAAAAAAAAAAAAAAA"}
    candidates = [
        scoring.evaluate(_primer_set(2, poor), profile),
        scoring.evaluate(_primer_set(1), profile),
    ]
    ranked = scoring.rank(candidates)
    assert [candidate.primer_set.index for candidate in ranked] == [1, 2]
    assert ranked[0].score.total < ranked[1].score.total


def test_rank_breaks_ties_by_the_designer_order():
    profile = profiles.load()
    candidates = [
        scoring.evaluate(_primer_set(3), profile),
        scoring.evaluate(_primer_set(1), profile),
    ]
    assert [item.primer_set.index for item in scoring.rank(candidates)] == [1, 3]
