import argparse
import math
import os
import random

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SEQ_MIN, SEQ_MAX = 1000, 3000
FRAG_MIN, FRAG_MAX = 100, 3000
N_SAMPLES = 10

RUN_MM = 80.0
SIZE_HIGH, SIZE_LOW = 4000, 50
LADDER = [100, 200, 300, 500, 750, 1000, 1500, 2000, 2500, 3000]


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


def take_samples(seq, n, rng):
    max_len = min(FRAG_MAX, len(seq))
    lengths = rng.sample(range(FRAG_MIN, max_len + 1), n)
    samples = []
    for length in lengths:
        start = rng.randint(0, len(seq) - length)
        samples.append((start, seq[start:start + length]))
    return samples


def migration_distance(length):
    length = min(max(length, SIZE_LOW), SIZE_HIGH)
    frac = (math.log10(SIZE_HIGH) - math.log10(length)) / (math.log10(SIZE_HIGH) - math.log10(SIZE_LOW))
    return RUN_MM * frac


def print_table(samples):
    rows = sorted(((migration_distance(len(f)), s, f) for s, f in samples))
    print(f"\n{'Band':<6}{'Start':>7}{'End':>7}{'Length (bp)':>13}{'Distance (mm)':>15}")
    print("-" * 48)
    for band, (dist, start, frag) in enumerate(rows, 1):
        print(f"{band:<6}{start + 1:>7}{start + len(frag):>7}{len(frag):>13}{dist:>15.1f}")
    print("(bands are numbered from the wells downward)")


def print_text_gel(samples, rows=40):
    lane = [""] * rows
    for _, frag in sorted(samples, key=lambda s: -len(s[1])):
        r = round(migration_distance(len(frag)) / RUN_MM * (rows - 1))
        while lane[r]:
            r += 1
        lane[r] = f"{len(frag):>4} bp"
    print("\n   wells")
    print("  +----------+")
    for row in lane:
        print(f"  | {row:<8} |")
    print("  +----------+")


def spread_labels(ys, min_gap):
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = list(ys)
    for a, b in zip(order, order[1:]):
        if out[b] - out[a] < min_gap:
            out[b] = out[a] + min_gap
    return out


def draw_gel(samples, title, outfile):
    fig, ax = plt.subplots(figsize=(5.5, 8))
    fig.patch.set_facecolor("#111111")
    ax.set_facecolor("#1d1d1d")

    ladder_x, sample_x, lane_w = 1.0, 3.0, 1.4
    ax.add_patch(Rectangle((0, -6), 4.4, RUN_MM + 12, color="#2b3a2b", zorder=0))

    for x in (ladder_x, sample_x):
        ax.add_patch(Rectangle((x - lane_w / 2, -5), lane_w, 4, color="#0a0a0a", zorder=1))

    for size in LADDER:
        y = migration_distance(size)
        ax.add_patch(Rectangle((ladder_x - lane_w / 2, y - 0.6), lane_w, 1.2, color="#9fd3ff", alpha=0.75, zorder=2))
        ax.text(ladder_x - lane_w / 2 - 0.1, y, f"{size}", color="#cccccc", fontsize=8, ha="right", va="center")

    ys = []
    for _, frag in samples:
        y = migration_distance(len(frag))
        ys.append(y)
        alpha = 0.35 + 0.65 * len(frag) / FRAG_MAX
        ax.add_patch(Rectangle((sample_x - lane_w / 2, y - 0.6), lane_w, 1.2, color="white", alpha=alpha, zorder=2))

    for y_text, y_band, (_, frag) in zip(spread_labels(ys, 2.6), ys, samples):
        ax.text(sample_x + lane_w / 2 + 0.1, y_text, f"{len(frag)} bp", color="#ffd479", fontsize=8, ha="left", va="center")

    ax.text(ladder_x, -7.5, "Ladder", color="white", ha="center", fontsize=10)
    ax.text(sample_x, -7.5, "Samples", color="white", ha="center", fontsize=10)
    ax.set_xlim(-0.6, 5.2)
    ax.set_ylim(RUN_MM + 8, -10)
    ax.axis("off")
    ax.set_title(title, color="white", fontsize=11, pad=10)
    fig.text(0.5, 0.02, "migration direction: top -> bottom (+ electrode at the bottom)", color="#aaaaaa", ha="center", fontsize=8)
    fig.savefig(outfile, dpi=150, facecolor=fig.get_facecolor())
    return fig


def main():
    parser = argparse.ArgumentParser(description="Gel electrophoresis simulation")
    parser.add_argument("fasta", nargs="?", default="sequence.fasta", help="FASTA file from NCBI (1000-3000 nt)")
    parser.add_argument("--seed", type=int, default=None, help="random seed for reproducible fragments")
    parser.add_argument("--start", type=int, default=0, help="start position if the sequence is longer than 3000 nt")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG but do not open a window")
    args = parser.parse_args()

    if not os.path.isfile(args.fasta):
        raise SystemExit(f"File not found: {args.fasta}")

    rng = random.Random(args.seed)
    header, seq = load_sequence(args.fasta, args.start)
    print(f"Sequence: {header[:70]}")
    print(f"Length used: {len(seq)} nt")

    samples = take_samples(seq, N_SAMPLES, rng)
    fragments = [frag for _, frag in samples]
    print(f"Took {len(fragments)} random fragments (lengths: {sorted(len(f) for f in fragments)})")

    print_table(samples)
    print_text_gel(samples)

    draw_gel(samples, f"Simulated agarose gel\n{header[:45]}", "ex1_gel.png")
    print("\nSaved gel image to ex1_gel.png")
    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()