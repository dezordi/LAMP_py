"""Tests for binding-site conservation across a target group."""

from __future__ import annotations

from lamppy import conservation
from lamppy.kmers import reverse_complement

SITE = "ACGTTGCATTCCGGATTACG"


def test_load_group_upper_cases_soft_masked_sequence(write_genome):
    genome = write_genome("g", {"c1": "acgtACGT"})
    group = conservation.load_group([genome])
    assert group == {"g.fasta": ("ACGTACGT",)}


def test_a_site_is_found_on_either_strand():
    assert conservation.holds_site([f"AAAA{SITE}TTTT"], SITE) is True
    assert (
        conservation.holds_site([f"AAAA{reverse_complement(SITE)}TTTT"], SITE) is True
    )
    assert conservation.holds_site(["AAAATTTT"], SITE) is False


def test_a_site_is_found_in_any_contig_of_a_genome():
    assert conservation.holds_site(["AAAA", f"TT{SITE}"], SITE) is True


def test_site_conservation_is_the_fraction_of_genomes_holding_the_site(write_genome):
    first = write_genome("g1", {"c1": f"AAAA{SITE}AAAA"})
    second = write_genome("g2", {"c1": f"CCCC{reverse_complement(SITE)}CCCC"})
    third = write_genome("g3", {"c1": "GGGGGGGGGGGG"})
    group = conservation.load_group([first, second, third])
    assert conservation.site_conservation(SITE, group) == 2 / 3


def test_soft_masking_never_hides_a_site(write_genome):
    genome = write_genome("g", {"c1": f"AAAA{SITE.lower()}AAAA"})
    group = conservation.load_group([genome])
    assert conservation.site_conservation(SITE, group) == 1.0


def test_an_empty_group_is_not_penalised():
    assert conservation.site_conservation(SITE, {}) == 1.0


def test_conservation_of_scores_every_oligo_of_a_set(write_genome):
    present, absent = SITE, "TTTTTTTTTTTTTTTTTTTT"
    genome = write_genome("g", {"c1": f"AAAA{present}AAAA"})
    group = conservation.load_group([genome])
    scored = conservation.conservation_of({"F3": present, "B3": absent}, group)
    assert scored == {"F3": 1.0, "B3": 0.0}
