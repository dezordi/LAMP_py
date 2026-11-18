"""Tests for the primer3 design backend."""

from __future__ import annotations

import random

import pytest

from lamppy import fasta, native, profiles, scoring
from lamppy.kmers import reverse_complement
from lamppy.primers import ORDER


@pytest.fixture
def profile() -> dict:
    return profiles.load()


@pytest.fixture
def interval() -> fasta.Fragment:
    generator = random.Random(7)
    sequence = "".join(generator.choice("ACGT") for _ in range(1500))
    return fasta.Fragment("chr1", 1000, sequence)


def test_the_oligo_of_a_forward_region_is_the_plus_strand():
    assert native.oligo("ACGTTTGGCC", 2, 4, "F3") == "GTTT"
    assert native.oligo("ACGTTTGGCC", 2, 4, "B1c") == "GTTT"


def test_the_oligo_of_a_reversed_region_faces_the_other_way():
    assert native.oligo("ACGTTTGGCC", 2, 4, "F1c") == reverse_complement("GTTT")
    assert native.oligo("ACGTTTGGCC", 2, 4, "B2") == reverse_complement("GTTT")
    assert native.oligo("ACGTTTGGCC", 2, 4, "B3") == reverse_complement("GTTT")


def test_the_oligo_is_upper_cased_so_soft_masking_never_reaches_a_primer():
    assert native.oligo("acgtttggcc", 2, 4, "F3") == "GTTT"


def test_sites_keep_only_temperatures_inside_the_window(
    interval: fasta.Fragment, profile: dict
):
    window = profiles.tm_range(profile, "F3")
    found = native.sites(interval.sequence, "F3", profile)
    assert found
    assert all(window[0] <= site.tm <= window[1] for site in found)
    assert all(17 <= site.length <= 22 for site in found)


def test_sites_keep_one_length_per_start(interval: fasta.Fragment, profile: dict):
    found = native.sites(interval.sequence, "F2", profile)
    starts = [site.start for site in found]
    assert len(starts) == len(set(starts))


def test_a_region_without_a_window_yields_no_sites(
    interval: fasta.Fragment, profile: dict
):
    assert native.sites(interval.sequence, "LF", profile) == []


def test_forward_arms_satisfy_the_configured_geometry(
    interval: fasta.Fragment, profile: dict
):
    geometry = profiles.parameters(profile)["geometry"]
    arms = list(native.forward_arms(interval.sequence, profile))
    assert arms
    for first, second, third in arms:
        assert geometry["f3_f2_gap"][0] <= second.start - first.end
        assert second.start - first.end <= geometry["f3_f2_gap"][1]
        assert geometry["f2_f1c"][0] <= third.end - second.start
        assert third.end - second.start <= geometry["f2_f1c"][1]
        assert third.start >= second.end


def test_backward_arms_satisfy_the_configured_geometry(
    interval: fasta.Fragment, profile: dict
):
    geometry = profiles.parameters(profile)["geometry"]
    arms = list(native.backward_arms(interval.sequence, profile))
    assert arms
    for first, second, third in arms:
        assert geometry["b1c_b2"][0] <= second.end - first.start
        assert second.end - first.start <= geometry["b1c_b2"][1]
        assert second.start >= first.end
        assert geometry["b2_b3_gap"][0] <= third.start - second.end
        assert third.start - second.end <= geometry["b2_b3_gap"][1]


def test_one_arm_is_kept_per_anchor_site(interval: fasta.Fragment, profile: dict):
    arms = list(native.forward_arms(interval.sequence, profile))
    anchors = [arm[1].start for arm in arms]
    assert len(anchors) == len(set(anchors))


def test_assembled_sets_respect_the_amplicon_window(
    interval: fasta.Fragment, profile: dict
):
    geometry = profiles.parameters(profile)["geometry"]
    assembled = native.assemble(interval.sequence, profile)
    assert assembled
    for chosen in assembled:
        span = chosen[-1].end - chosen[0].start
        assert geometry["amplicon"][0] <= span <= geometry["amplicon"][1]
        assert [site.start for site in chosen] == sorted(site.start for site in chosen)


def test_assembly_stops_at_the_limit(interval: fasta.Fragment, profile: dict):
    assert len(native.assemble(interval.sequence, profile, limit=5)) == 5


def test_a_sequence_too_short_for_a_set_yields_nothing(profile: dict):
    assert native.assemble("ACGT" * 20, profile) == []


def test_tm_deviation_rewards_temperatures_near_the_middle_of_the_window(
    interval: fasta.Fragment, profile: dict
):
    assembled = native.assemble(interval.sequence, profile)
    deviations = [native.tm_deviation(chosen, profile) for chosen in assembled]
    assert min(deviations) >= 0.0
    assert min(deviations) < max(deviations)


def test_design_fragment_reports_sets_in_source_coordinates(
    interval: fasta.Fragment, profile: dict
):
    designed = native.design_fragment(interval, profile, sets=2)
    assert designed
    first = designed[0]
    assert first.reference == "chr1:1000"
    assert [primer.name for primer in first.primers] == list(ORDER)
    assert all(primer.contig == "chr1" for primer in first.primers)
    assert all(primer.start >= 1000 for primer in first.primers)


def test_a_designed_set_sits_where_its_sequence_really_is(
    interval: fasta.Fragment, profile: dict
):
    first = native.design_fragment(interval, profile, sets=1)[0]
    for primer in first.primers:
        offset = primer.start - interval.start
        window = interval.sequence[offset : offset + len(primer)].upper()
        assert primer.sequence in (window, reverse_complement(window))


def test_design_fragment_over_generates_for_the_ranking(
    interval: fasta.Fragment, profile: dict
):
    designed = native.design_fragment(interval, profile, sets=3)
    assert len(designed) > 3
    assert len({primer_set.index for primer_set in designed}) == len(designed)


def test_a_designed_set_can_be_scored_like_any_other(
    interval: fasta.Fragment, profile: dict
):
    designed = native.design_fragment(interval, profile, sets=2)
    ranked = scoring.rank(scoring.evaluate(item, profile) for item in designed)
    assert ranked[0].score.total <= ranked[-1].score.total
    assert set(ranked[0].primer_set.reaction_oligos()) == {"F3", "B3", "FIP", "BIP"}


def test_designing_without_a_profile_uses_the_default(interval: fasta.Fragment):
    assert native.design_fragment(interval, sets=1)
