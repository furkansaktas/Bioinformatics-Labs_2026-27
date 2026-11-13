import argparse
import math
import os
from collections import Counter
from itertools import combinations

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ENZYMES = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI": ("TCGA", 1),
    "HaeIII": ("GGCC", 2),
}

RUN_MM = 80.0
SIZE_LOW = 10
PAIR = dict(zip("ACGT", "TGCA"))


def revcomp(seq):
    return "".join(PAIR[b] for b in reversed(seq))


def read_first_fasta_record(path):
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            elif header is not None:
                chunks.append(line.upper())
    if header is None:
        raise SystemExit("No FASTA record found (the file should start with '>').")
    return header, "".join(chunks)


def find_all(text, pattern):
    positions, i = [], text.find(pattern)
    while i != -1:
        positions.append(i)
        i = text.find(pattern, i + 1)
    return positions


def find_cuts(seq, site, offset, circular):
    length = len(seq)
    text = seq + seq[:len(site) - 1] if circular else seq
    cuts = [p + offset for p in find_all(text, site) if p < length]
    rc = revcomp(site)
    if rc != site:
        cuts += [p + len(site) - offset for p in find_all(text, rc) if p < length]
    if circular:
        cuts = [c % length for c in cuts]
    return sorted(set(cuts))


def make_fragments(cuts, length, circular):
    cuts = sorted(set(cuts))
    if not cuts:
        return [(1, length, length)]
    if circular:
        pieces = [(cuts[i] + 1, cuts[i + 1], cuts[i + 1] - cuts[i]) for i in range(len(cuts) - 1)]
        pieces.append((cuts[-1] + 1, cuts[0], length - cuts[-1] + cuts[0]))
        return pieces
    bounds = [0] + cuts + [length]
    return [(bounds[i] + 1, bounds[i + 1], bounds[i + 1] - bounds[i]) for i in range(len(bounds) - 1)]


def combined_cuts(names, cuts_by_enzyme):
    return sorted(set(c for n in names for c in cuts_by_enzyme[n]))


def cut_marker(site, offset):
    return f"{site[:offset]}^{site[offset:]}"


def print_enzyme_table(enzymes, cuts_by_enzyme):
    cut_header = "Cut (5'->3')"
    print(f"\n{'Enzyme':<10}{'Recognized site':<18}{cut_header:<16}{'Cleavages':>10}")
    print("-" * 54)
    for name, (site, offset) in enzymes.items():
        print(f"{name:<10}{site:<18}{cut_marker(site, offset):<16}{len(cuts_by_enzyme[name]):>10}")


def print_combinations(enzymes, cuts_by_enzyme, length, circular):
    names = list(enzymes)
    print(f"\n{'Combination':<44}{'Cuts':>6}{'Fragments':>11}{'Smallest':>10}{'Largest':>10}")
    print("-" * 81)
    for r in range(1, len(names) + 1):
        for combo in combinations(names, r):
            cuts = combined_cuts(combo, cuts_by_enzyme)
            sizes = [f[2] for f in make_fragments(cuts, length, circular)]
            print(f"{' + '.join(combo):<44}{len(cuts):>6}{len(sizes):>11}{min(sizes):>10,}{max(sizes):>10,}")


def shorten(items, limit, full):
    if full or len(items) <= limit:
        return items, False
    return items[:limit - 3] + items[-3:], True


def print_detail(names, cuts_by_enzyme, length, circular, full):
    cuts = combined_cuts(names, cuts_by_enzyme)
    fragments = make_fragments(cuts, length, circular)
    print(f"\nDigest with: {' + '.join(names)}")
    print("=" * 60)
    print("Cleavage positions (the cut is made after this base):")
    for n in names:
        shown = ", ".join(f"{c:,}" for c in cuts_by_enzyme[n][:12])
        more = " ..." if len(cuts_by_enzyme[n]) > 12 and not full else ""
        if full:
            shown = ", ".join(f"{c:,}" for c in cuts_by_enzyme[n])
        print(f"  {n} ({len(cuts_by_enzyme[n])}): {shown if shown else '-'}{more}")
    print(f"\nTotal cleavages: {len(cuts)}   Fragments: {len(fragments)}")
    print(f"\n{'#':<6}{'Start':>10}{'End':>10}{'Length (bp)':>14}")
    print("-" * 40)
    numbered = list(enumerate(fragments, 1))
    rows, cut_short = shorten(numbered, 20, full)
    for idx, (start, end, size) in rows:
        print(f"{idx:<6}{start:>10,}{end:>10,}{size:>14,}")
        if cut_short and idx == 17:
            print(f"{'...':<6}{'(use --full to list all fragments)':>34}")
    print("-" * 40)
    print(f"{'Sum':<6}{'':>20}{sum(f[2] for f in fragments):>14,}")


def migration_distance(size, size_high):
    size = min(max(size, SIZE_LOW), size_high)
    frac = (math.log10(size_high) - math.log10(size)) / (math.log10(size_high) - math.log10(SIZE_LOW))
    return RUN_MM * frac


def ladder_sizes(size_high):
    return [m * 10 ** k for k in range(2, 9) for m in (1, 2, 5) if SIZE_LOW < m * 10 ** k < size_high]


def size_label(size):
    return f"{size / 1000:g} kb" if size >= 1000 else f"{size} bp"


def draw_gel(lanes, ladder, length, title, outfile):
    size_high = length * 1.3
    fig, ax = plt.subplots(figsize=(1.4 * (len(lanes) + 2) + 1.5, 8))
    fig.patch.set_facecolor("#111111")
    ax.set_facecolor("#1d1d1d")
    spacing, lane_w = 1.5, 1.1
    total_w = spacing * (len(lanes) + 2)
    ax.add_patch(Rectangle((0, -6), total_w, RUN_MM + 12, color="#2b3a2b", zorder=0))

    all_lanes = [("Ladder", None)] + lanes
    for k, (name, sizes) in enumerate(all_lanes):
        x = spacing * (k + 1)
        ax.add_patch(Rectangle((x - lane_w / 2, -5), lane_w, 4, color="#0a0a0a", zorder=1))
        ax.text(x, -7.5, name, color="white", ha="center", va="bottom", fontsize=9)
        if sizes is None:
            for size in ladder:
                y = migration_distance(size, size_high)
                ax.add_patch(Rectangle((x - lane_w / 2, y - 0.5), lane_w, 1.0, color="#9fd3ff", alpha=0.8, zorder=2))
                ax.text(x - lane_w / 2 - 0.08, y, size_label(size), color="#cccccc", fontsize=7, ha="right", va="center")
            continue
        counts = Counter(sizes)
        top_mass = max(size * n for size, n in counts.items())
        for size, n in counts.items():
            alpha = 0.2 + 0.8 * math.sqrt(size * n / top_mass)
            y = migration_distance(size, size_high)
            ax.add_patch(Rectangle((x - lane_w / 2, y - 0.5), lane_w, 1.0, color="white", alpha=alpha, zorder=2))

    ax.set_xlim(-0.6, total_w + 0.3)
    ax.set_ylim(RUN_MM + 8, -12)
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=11, pad=10)
    fig.text(0.5, 0.015, "migration direction: top -> bottom (fragments under 10 bp are drawn at the bottom edge)",
             color="#aaaaaa", ha="center", fontsize=8)
    fig.savefig(outfile, dpi=150, facecolor=fig.get_facecolor())
    return fig


def main():
    parser = argparse.ArgumentParser(description="Restriction enzyme digest simulator")
    parser.add_argument("fasta", nargs="?", default="sequence.fasta", help="FASTA file (bacterial or viral genome)")
    parser.add_argument("--enzymes", nargs="+", default=None, help="enzymes to combine in the detailed report and gel")
    parser.add_argument("--add-enzyme", nargs=3, action="append", metavar=("NAME", "SITE", "CUT"),
                        help="extra enzyme: name, recognized sequence, cut position inside the site (e.g. SmaI CCCGGG 3)")
    parser.add_argument("--circular", action="store_true", help="treat the DNA as circular (bacterial chromosome or plasmid)")
    parser.add_argument("--full", action="store_true", help="print every cleavage position and fragment")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG but do not open a window")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta):
        raise SystemExit(f"File not found: {args.fasta}")

    enzymes = dict(ENZYMES)
    for name, site, cut in args.add_enzyme or []:
        site = site.upper()
        if not set(site) <= set("ACGT") or not cut.isdigit() or not 1 <= int(cut) < len(site):
            raise SystemExit(f"Invalid enzyme definition: {name} {site} {cut}")
        enzymes[name] = (site, int(cut))

    selected = args.enzymes or list(enzymes)
    unknown = [n for n in selected if n not in enzymes]
    if unknown:
        raise SystemExit(f"Unknown enzyme(s): {', '.join(unknown)}. Available: {', '.join(enzymes)}")

    header, seq = read_first_fasta_record(args.fasta)
    length = len(seq)
    print(f"Sequence: {header[:70]}")
    print(f"Length: {length:,} bp ({'circular' if args.circular else 'linear'})")

    cuts_by_enzyme = {n: find_cuts(seq, site, offset, args.circular) for n, (site, offset) in enzymes.items()}

    print_enzyme_table(enzymes, cuts_by_enzyme)
    print_combinations(enzymes, cuts_by_enzyme, length, args.circular)
    print_detail(selected, cuts_by_enzyme, length, args.circular, args.full)

    lanes = [(n, [f[2] for f in make_fragments(cuts_by_enzyme[n], length, args.circular)]) for n in enzymes]
    combo_cuts = combined_cuts(selected, cuts_by_enzyme)
    lanes.append(("Combined", [f[2] for f in make_fragments(combo_cuts, length, args.circular)]))
    title = f"Simulated agarose gel\n{header[:48]}\nCombined = {' + '.join(selected)}"
    draw_gel(lanes, ladder_sizes(length * 1.3), length, title, "lab7_gel.png")
    print("\nSaved gel image to lab7_gel.png")
    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()