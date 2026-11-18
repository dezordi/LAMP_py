"""Tests for oligo and set thermodynamics."""

from __future__ import annotations

from lamppy import thermo

F3 = "CGGTTCCCTGTACTCGAA"
F2 = "CGAAATCTTCAAACACTACG"
F1C = "TGCTAATTCCTTTGTTGAGGCCG"
B1C = "CCTGACGGAAGCAGCATTAAG"
B2 = "TCGAACGTAACCAAAGTCGT"
B3 = "ACGTTGGCATTCGATTCC"
REGIONS = {"F3": F3, "F2": F2, "F1c": F1C, "B1c": B1C, "B2": B2, "B3": B3}


def test_gc_fraction_counts_both_cases():
    assert thermo.gc_fraction("acGT") == 0.5
    assert thermo.gc_fraction("") == 0.0


def test_melting_temperature_rises_with_gc_content():
    assert thermo.melting_temperature(
        "ATATATATATATATATAT"
    ) < thermo.melting_temperature("GCGCGCGCGCGCGCGCGC")


def test_melting_temperature_ignores_case():
    assert thermo.melting_temperature(F2.lower()) == thermo.melting_temperature(F2)


def test_hairpin_is_negative_for_a_self_complementary_oligo():
    assert thermo.hairpin_dg("GGGGGGCCCCTTTTGGGGCCCCCC") < -5.0


def test_dimer_is_strongly_negative_for_complementary_oligos():
    assert thermo.dimer_dg("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT") < -10.0


def test_end_stability_measures_only_the_terminal_bases():
    assert -6.0 < thermo.end_stability(F2) < -1.0


def test_end_stability_of_the_five_prime_end_differs_from_the_three_prime_end():
    assert thermo.end_stability(F1C, five_prime=True) != thermo.end_stability(F1C)


def test_end_stability_of_a_gc_end_is_more_stable_than_an_at_end():
    assert thermo.end_stability("TTTTTTTTTTTTTTTGCGCG") < thermo.end_stability(
        "GCGCGCGCGCGCGCGATATA"
    )


def test_evaluate_primer_reports_every_measurement():
    measured = thermo.evaluate_primer("F2", F2)
    assert measured.name == "F2"
    assert 40.0 < measured.tm < 70.0
    assert measured.gc == thermo.gc_fraction(F2)
    assert measured.end_stability < 0.0


def test_evaluate_set_describes_every_region():
    measured = thermo.evaluate_set(REGIONS)
    assert [primer.name for primer in measured.primers] == list(REGIONS)
    assert measured.primer("F1c") is not None
    assert measured.primer("LF") is None


def test_evaluate_set_reports_the_offset_between_the_inner_regions():
    measured = thermo.evaluate_set(REGIONS)
    expected = thermo.melting_temperature(F1C) - thermo.melting_temperature(F2)
    assert measured.offsets["F2_F1c"] == expected
    assert set(measured.offsets) == {"F2_F1c", "B2_B1c"}


def test_evaluate_set_measures_duplexes_over_the_oligos_given():
    oligos = {"F3": F3, "B3": B3, "FIP": F1C + F2, "BIP": B1C + B2}
    composite = thermo.evaluate_set(REGIONS, oligos)
    regions_only = thermo.evaluate_set(REGIONS)
    assert composite.worst_dimer != regions_only.worst_dimer
    assert composite.primers == regions_only.primers


def test_the_worst_end_dimer_is_not_more_stable_than_the_worst_dimer():
    measured = thermo.evaluate_set(REGIONS)
    assert measured.worst_dimer <= measured.worst_end_dimer


def test_an_empty_set_has_no_interactions():
    measured = thermo.evaluate_set({})
    assert measured.primers == ()
    assert measured.worst_dimer == 0.0
    assert measured.worst_end_dimer == 0.0
    assert measured.worst_hairpin == 0.0
