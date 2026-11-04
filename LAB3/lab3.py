import argparse
import os
from collections import Counter

import matplotlib.pyplot as plt

DEFAULT_FASTA = "GCF_009858895.2_ASM985889v3_genomic.fna"

BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: AMINO_ACIDS[i]
    for i, (a, b, c) in enumerate((a, b, c) for a in BASES for b in BASES for c in BASES)
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    """Return {header: sequence} for every record in a FASTA file."""
    records, header, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line.upper())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def count_codons_six_frames(seq):
    """Count every complete, unambiguous codon in all 6 reading frames."""
    counts = Counter()
    for strand in (seq, reverse_complement(seq)):
        for frame in range(3):
            for i in range(frame, len(strand) - 2, 3):
                codon = strand[i:i + 3]
                if codon in CODON_TABLE:        # skips codons containing N, etc.
                    counts[codon] += 1
    return counts


def amino_acid_usage(codon_counts):
    """Collapse codon counts into amino-acid counts."""
    usage = Counter()
    for codon, n in codon_counts.items():
        usage[CODON_TABLE[codon]] += n
    return usage


def print_report(header, counts, top_n):
    total = sum(counts.values())
    print("\n" + "=" * 60)
    print(f" {header[:56]}")
    print(f" Total codons (6 frames): {total:,}")
    print("=" * 60)
    print(f" {'Rank':<5}{'Codon':<8}{'AA':<5}{'Count':>10}{'Freq %':>10}")
    print("-" * 60)
    for rank, (codon, n) in enumerate(counts.most_common(top_n), 1):
        print(f" {rank:<5}{codon:<8}{CODON_TABLE[codon]:<5}{n:>10,}{n / total * 100:>10.3f}")

    usage = amino_acid_usage(counts)
    print("-" * 60)
    print(" Top 5 amino acids (all synonymous codons combined):")
    for aa, n in usage.most_common(5):
        print(f"   {aa}: {n:,} ({n / total * 100:.2f}%)")
    stops = usage.get("*", 0)
    print(f" Stop codons: {stops:,} ({stops / total * 100:.2f}%)")
    print("=" * 60)


def plot_top_codons(header, counts, top_n, outfile):
    total = sum(counts.values())
    top = counts.most_common(top_n)
    labels = [f"{c}\n{CODON_TABLE[c]}" for c, _ in top]
    freqs = [n / total * 100 for _, n in top]

    # same colour for synonymous codons
    aas = [CODON_TABLE[c] for c, _ in top]
    palette = plt.colormaps["tab20"]
    colour_of = {aa: palette(i % 20) for i, aa in enumerate(dict.fromkeys(aas))}

    fig, ax = plt.subplots(figsize=(11, 5))
    bars = ax.bar(range(len(top)), freqs, color=[colour_of[a] for a in aas], edgecolor="black")
    for bar, f in zip(bars, freqs):
        ax.text(bar.get_x() + bar.get_width() / 2, f + 0.02, f"{f:.2f}%", ha="center", fontsize=8)

    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(labels)
    ax.set_xlabel("Codon (amino acid)")
    ax.set_ylabel("Frequency (%)")
    ax.set_title(f"Top {top_n} codons, 6 reading frames\n{header[:70]}")
    ax.set_ylim(0, max(freqs) * 1.15)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Six-frame codon frequency analyzer")
    parser.add_argument("fasta", nargs="?", default=DEFAULT_FASTA, help="FASTA file path")
    parser.add_argument("-n", "--top", type=int, default=10, help="how many codons to show")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG but don't open a window")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta):
        raise SystemExit(f"File not found: {args.fasta}")

    records = read_fasta(args.fasta)
    if not records:
        raise SystemExit("No sequences found in the FASTA file.")
    print(f"Read {len(records)} sequence(s) from {args.fasta}")

    for i, (header, seq) in enumerate(records.items(), 1):
        print(f"\nProcessing record {i}: {len(seq):,} bp")
        counts = count_codons_six_frames(seq)
        print_report(header, counts, args.top)
        plot_top_codons(header, counts, args.top, f"codon_plot_{i}.png")
        print(f"Saved plot to codon_plot_{i}.png")

    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()