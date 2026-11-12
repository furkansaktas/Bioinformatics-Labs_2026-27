import argparse
import os
from collections import Counter

import matplotlib.pyplot as plt

SEQ_MIN, SEQ_MAX = 1000, 3000
UNIT_MIN, UNIT_MAX = 3, 5
REPEAT_MIN, REPEAT_MAX = 6, 10
VALID = set("ACGT")


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


def load_sequence(path, start):
    header, seq = read_first_fasta_record(path)
    if len(seq) < SEQ_MIN:
        raise SystemExit(f"Sequence is only {len(seq)} nt; the exercise needs {SEQ_MIN}-{SEQ_MAX} nt.")
    if len(seq) > SEQ_MAX:
        print(f"Note: sequence has {len(seq):,} nt, using {SEQ_MAX} nt starting at position {start}.")
        seq = seq[start:start + SEQ_MAX]
        if len(seq) < SEQ_MIN:
            raise SystemExit("--start leaves fewer than 1000 nt; choose a smaller value.")
    return header, seq


def is_primitive(unit):
    n = len(unit)
    return not any(n % d == 0 and unit[:d] * (n // d) == unit for d in range(1, n))


def find_tandem_repeats(seq):
    repeats, too_long = [], 0
    for k in range(UNIT_MIN, UNIT_MAX + 1):
        i = 0
        while i <= len(seq) - 2 * k:
            unit = seq[i:i + k]
            if not set(unit) <= VALID or not is_primitive(unit):
                i += 1
                continue

            copies = 1
            while seq[i + copies * k:i + (copies + 1) * k] == unit:
                copies += 1

            if copies >= 2:
                length = copies * k
                if REPEAT_MIN <= length <= REPEAT_MAX:
                    repeats.append({"start": i + 1, "end": i + length, "unit": unit,
                                    "copies": copies, "length": length})
                elif length > REPEAT_MAX:
                    too_long += 1
                i += length
            else:
                i += 1
    repeats.sort(key=lambda r: (r["start"], r["unit"]))
    return repeats, too_long


def print_report(header, seq, repeats, too_long):
    print(f"\nSequence: {header[:70]}")
    print(f"Length: {len(seq)} nt")
    print(f"Searching for repeats: unit {UNIT_MIN}-{UNIT_MAX} bases, region {REPEAT_MIN}-{REPEAT_MAX} bases\n")

    print(f"{'#':<4}{'Start':>6}{'End':>7}  {'Unit':<7}{'Copies':>7}{'Length':>8}  Repeat")
    print("-" * 62)
    for n, r in enumerate(repeats, 1):
        print(f"{n:<4}{r['start']:>6}{r['end']:>7}  {r['unit']:<7}{r['copies']:>7}{r['length']:>8}  {r['unit'] * r['copies']}")

    print("-" * 62)
    print(f"Total repeats found: {len(repeats)}")
    for k in range(UNIT_MIN, UNIT_MAX + 1):
        print(f"  unit length {k}: {sum(1 for r in repeats if len(r['unit']) == k)}")
    common = Counter(r["unit"] for r in repeats).most_common(3)
    if common:
        print("Most common units: " + ", ".join(f"{u} (x{c})" for u, c in common))
    if too_long:
        print(f"Skipped {too_long} longer repeat region(s) (> {REPEAT_MAX} bases)")


def plot_repeats(header, seq, repeats, outfile):
    colours = {3: "#1f77b4", 4: "#ff7f0e", 5: "#2ca02c"}
    fig, ax = plt.subplots(figsize=(11, 3.5))
    for k in range(UNIT_MIN, UNIT_MAX + 1):
        pts = [r for r in repeats if len(r["unit"]) == k]
        ax.scatter([r["start"] for r in pts], [k] * len(pts), s=40, color=colours[k],
                   label=f"{k}-base unit ({len(pts)})", zorder=3)
    ax.set_xlim(0, len(seq))
    ax.set_ylim(UNIT_MIN - 0.7, UNIT_MAX + 0.7)
    ax.set_yticks(range(UNIT_MIN, UNIT_MAX + 1))
    ax.set_xlabel("Position in sequence (nt)")
    ax.set_ylabel("Repeat unit length")
    ax.set_title(f"Tandem repeats ({REPEAT_MIN}-{REPEAT_MAX} bases)\n{header[:70]}")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Tandem repeat finder")
    parser.add_argument("fasta", nargs="?", default="sequence.fasta", help="FASTA file from NCBI (1000-3000 nt)")
    parser.add_argument("--start", type=int, default=0, help="start position if the sequence is longer than 3000 nt")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG but do not open a window")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta):
        raise SystemExit(f"File not found: {args.fasta}")

    header, seq = load_sequence(args.fasta, args.start)
    repeats, too_long = find_tandem_repeats(seq)
    print_report(header, seq, repeats, too_long)

    plot_repeats(header, seq, repeats, "ex2_repeats.png")
    print("\nSaved plot to ex2_repeats.png")
    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()