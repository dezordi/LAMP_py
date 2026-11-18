"""Canonical k-mer counting and cross-genome sharing statistics.

Counts are held in memory, which is adequate for bacterial and viral genome
sets. KMC3 or unikmer would replace this module if eukaryotic inputs are ever
in scope.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Iterator, Mapping
from pathlib import Path

from lamppy import fasta

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_UNKNOWN = 0


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an ACGT sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its complement.

    Collapsing both strands into one key makes a count strand-agnostic, which
    is what primer binding requires: a site present on the reverse strand of a
    background genome is just as disqualifying as one on the forward strand.
    """
    upper = kmer.upper()
    reverse = reverse_complement(upper)
    return min(upper, reverse)


def iter_kmers(sequence: str, size: int) -> Iterator[tuple[int, str]]:
    """Yield (offset, canonical k-mer) for every window of the sequence."""
    for start in range(len(sequence) - size + 1):
        yield start, canonical(sequence[start : start + size])


def count_genome(path: Path, size: int) -> Counter[str]:
    """Count canonical k-mers across the unambiguous fragments of a genome."""
    counts: Counter[str] = Counter()
    for fragment in fasta.iter_fragments(path, min_length=size):
        counts.update(kmer for _, kmer in iter_kmers(fragment.sequence, size))
    return counts


def count_genomes(paths: Iterable[Path], size: int) -> dict[str, Counter[str]]:
    """Count canonical k-mers per genome, keyed by file name.

    File names rather than stems are the keys, so two inputs differing only in
    extension remain distinct genomes.
    """
    return {path.name: count_genome(path, size) for path in paths}


def genome_frequency(counts: Mapping[str, Counter[str]]) -> Counter[str]:
    """Count in how many genomes each canonical k-mer occurs.

    Occurrences inside one genome collapse to one, so the value is a genome
    count and never a copy number.
    """
    frequency: Counter[str] = Counter()
    for genome_counts in counts.values():
        frequency.update(genome_counts.keys())
    return frequency


def sharing_rows(
    counts: Mapping[str, Counter[str]], min_genomes: int = 2
) -> Iterator[tuple[str, int, int, tuple[int, ...]]]:
    """Yield (k-mer, genomes sharing it, total occurrences, per-genome counts).

    Rows are restricted to k-mers seen in at least min_genomes genomes and are
    ordered by decreasing sharing, then decreasing total occurrences.
    """
    labels = list(counts)
    totals = genome_frequency(counts)
    shared = (kmer for kmer, genomes in totals.items() if genomes >= min_genomes)
    rows = (
        (
            kmer,
            totals[kmer],
            sum(counts[label].get(kmer, _UNKNOWN) for label in labels),
            tuple(counts[label].get(kmer, _UNKNOWN) for label in labels),
        )
        for kmer in shared
    )
    yield from sorted(rows, key=lambda row: (-row[1], -row[2], row[0]))


def write_sharing_table(
    counts: Mapping[str, Counter[str]], path: Path, min_genomes: int = 2
) -> int:
    """Write the k-mer sharing table as TSV and return how many rows it holds."""
    labels = list(counts)
    written = 0
    with path.open("w") as handle:
        handle.write("kmer\tn_genomes\ttotal_count\t" + "\t".join(labels) + "\n")
        for kmer, genomes, total, per_genome in sharing_rows(counts, min_genomes):
            values = "\t".join(str(value) for value in per_genome)
            handle.write(f"{kmer}\t{genomes}\t{total}\t{values}\n")
            written += 1
    return written


def spectrum(counts: Mapping[str, Counter[str]]) -> dict[int, int]:
    """Count distinct k-mers by how many genomes hold them, keyed by that number.

    The shape of this distribution is the cheapest feasibility check on a genome
    set: most k-mers held by every genome means a near-clonal set with no private
    sequence to target, while a large share held by one genome means there is
    signal to work with.
    """
    frequency = genome_frequency(counts)
    return dict(sorted(Counter(frequency.values()).items()))


def write_spectrum(counts: Mapping[str, Counter[str]], path: Path) -> int:
    """Write the k-mer spectrum as TSV and return how many rows it holds."""
    distribution = spectrum(counts)
    total = sum(distribution.values())
    with path.open("w") as handle:
        handle.write("n_genomes\tn_kmers\tfraction\n")
        for genomes, kmers_held in distribution.items():
            handle.write(f"{genomes}\t{kmers_held}\t{kmers_held / total:.6f}\n")
    return len(distribution)


def sharing_matrix(
    counts: Mapping[str, Counter[str]],
) -> list[tuple[str, str, int, float]]:
    """Return (genome, genome, shared distinct k-mers, Jaccard) for every pair.

    The Jaccard index over k-mer sets is the neighbour distance that decides
    which genomes belong in a background set, so no separate sketching tool is
    needed to obtain it.
    """
    labels = list(counts)
    sets = {label: set(counts[label]) for label in labels}
    return [
        (
            first,
            second,
            len(sets[first] & sets[second]),
            (
                len(sets[first] & sets[second]) / len(sets[first] | sets[second])
                if sets[first] or sets[second]
                else 0.0
            ),
        )
        for index, first in enumerate(labels)
        for second in labels[index + 1 :]
    ]


def write_sharing_matrix(counts: Mapping[str, Counter[str]], path: Path) -> int:
    """Write pairwise k-mer sharing as TSV and return how many pairs it holds.

    One row per unordered pair rather than a square matrix, so the table feeds a
    heatmap or a clustering directly without reshaping.
    """
    pairs = sharing_matrix(counts)
    with path.open("w") as handle:
        handle.write("genome_a\tgenome_b\tshared_kmers\tjaccard\n")
        for first, second, shared, jaccard in pairs:
            handle.write(f"{first}\t{second}\t{shared}\t{jaccard:.6f}\n")
    return len(pairs)


def summary_rows(
    counts: Mapping[str, Counter[str]],
) -> Iterator[tuple[str, int, int, int, int, int, float]]:
    """Yield one per-genome summary of k-mer totals, sharing and copy number.

    The fields are the genome, its total and distinct k-mers, how many are
    private to it, how many it shares with another genome, its highest
    intra-genome copy number, and the fraction of its k-mer positions whose
    k-mer occurs more than once in it.
    """
    frequency = genome_frequency(counts)
    for label, genome_counts in counts.items():
        total = sum(genome_counts.values())
        private = sum(1 for kmer in genome_counts if frequency[kmer] == 1)
        multicopy = sum(value for value in genome_counts.values() if value > 1)
        yield (
            label,
            total,
            len(genome_counts),
            private,
            len(genome_counts) - private,
            max(genome_counts.values(), default=0),
            multicopy / total if total else 0.0,
        )


def write_summary(counts: Mapping[str, Counter[str]], path: Path) -> int:
    """Write the per-genome k-mer summary as TSV and return how many rows it holds."""
    written = 0
    with path.open("w") as handle:
        handle.write(
            "genome\ttotal_kmers\tdistinct_kmers\tprivate_kmers\tshared_kmers\t"
            "max_self_copies\tmulticopy_fraction\n"
        )
        for row in summary_rows(counts):
            handle.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\t"
                f"{row[6]:.6f}\n"
            )
            written += 1
    return written
