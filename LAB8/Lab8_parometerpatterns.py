import argparse
import csv
import os

import matplotlib.pyplot as plt

BOARD_SEQUENCE = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTAT" + "CTAGCGATCTA"
WINDOW = 30
STEP = 1
MAX_PATTERN_PLOTS = 4
PREVIEW_ROWS = 20


def read_fasta(path):
    records, header, chunks = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line.upper())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def gc_content(window):
    return 100 * sum(base in "GC" for base in window) / len(window)


def information_content(window):
    n = len(window)
    total = 0.0
    for lag in range(1, n):
        matches = sum(window[i] == window[i + lag] for i in range(n - lag))
        total += matches / (n - lag)
    return 100 - 100 * total / (n - 1)


def sliding_pattern(seq, window=WINDOW, step=STEP):
    pattern = []
    for start in range(0, len(seq) - window + 1, step):
        piece = seq[start:start + window]
        pattern.append({
            "start": start + 1,
            "end": start + window,
            "center": start + (window + 1) / 2,
            "gc": gc_content(piece),
            "info": information_content(piece),
        })
    return pattern


def center_of_weight(positions, values):
    total = sum(values)
    if total == 0:
        return sum(positions) / len(positions)
    return sum(p * v for p, v in zip(positions, values)) / total


def analyse(name, seq, window, step):
    pattern = sliding_pattern(seq, window, step)
    positions = [w["center"] for w in pattern]
    gc_values = [w["gc"] for w in pattern]
    info_values = [w["info"] for w in pattern]
    return {
        "name": name,
        "seq": seq,
        "pattern": pattern,
        "positions": positions,
        "gc_values": gc_values,
        "info_values": info_values,
        "gc_center": center_of_weight(positions, gc_values),
        "info_center": center_of_weight(positions, info_values),
    }


def print_window_table(result, full):
    pattern = result["pattern"]
    print(f"\n{'Window':<8}{'Start':>7}{'End':>7}{'Center':>9}{'C+G %':>9}{'Info':>9}")
    print("-" * 49)
    rows = list(enumerate(pattern, 1))
    if not full and len(rows) > PREVIEW_ROWS:
        rows = rows[:PREVIEW_ROWS - 4] + [None] + rows[-3:]
    for row in rows:
        if row is None:
            print(f"{'...':<8}{'(use --full to list every window)':>41}")
            continue
        idx, w = row
        print(f"{idx:<8}{w['start']:>7}{w['end']:>7}{w['center']:>9.1f}{w['gc']:>9.2f}{w['info']:>9.2f}")


def print_summary(result, window):
    seq, name = result["seq"], result["name"]
    n_windows = len(result["pattern"])
    mean_gc = sum(result["gc_values"]) / n_windows
    mean_info = sum(result["info_values"]) / n_windows
    print(f"\nSequence: {name[:70]}")
    print(f"Length: {len(seq)} bp, window {window} bp, {n_windows} windows")
    print(f"Whole sequence:  C+G = {gc_content(seq):.2f}%   information content = {information_content(seq):.2f}")
    print(f"Window average:  C+G = {mean_gc:.2f}%   information content = {mean_info:.2f}")
    print(f"Window maximum:  C+G = {max(result['gc_values']):.2f}%   information content = {max(result['info_values']):.2f}")
    print(f"Window minimum:  C+G = {min(result['gc_values']):.2f}%   information content = {min(result['info_values']):.2f}")
    for label, key in (("C+G pattern", "gc_center"), ("Information pattern", "info_center")):
        c = result[key]
        print(f"Center of weight, {label}: position {c:.2f} bp ({100 * c / len(seq):.1f}% of the sequence)")


def draw_pattern(result, outfile):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
    x = result["positions"]
    ax1.plot(x, result["gc_values"], color="#1f77b4", linewidth=2)
    ax1.axvline(result["gc_center"], color="#1f77b4", linestyle="--", alpha=0.7,
                label=f"center of weight = {result['gc_center']:.1f}")
    ax1.set_ylabel("C+G (%)")
    ax1.legend(loc="upper left", fontsize=9)
    ax1.grid(True, linestyle="--", alpha=0.4)

    ax2.plot(x, result["info_values"], color="#ff7f0e", linewidth=2)
    ax2.axvline(result["info_center"], color="#ff7f0e", linestyle="--", alpha=0.7,
                label=f"center of weight = {result['info_center']:.1f}")
    ax2.set_ylabel("Information content")
    ax2.set_xlabel("Window center position (bp)")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.grid(True, linestyle="--", alpha=0.4)

    ax1.set_title(f"Pattern of {result['name'][:60]}")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    return fig


def draw_centers(results, outfile):
    fig, ax = plt.subplots(figsize=(10, max(2.5, 0.6 * len(results) + 1.5)))
    for row, r in enumerate(results):
        lo, hi = r["positions"][0], r["positions"][-1]
        ax.hlines(row, lo, hi, color="#cccccc", linewidth=6, zorder=1)
        ax.scatter(r["gc_center"], row, color="#1f77b4", s=90, zorder=3, label="C+G pattern" if row == 0 else None)
        ax.scatter(r["info_center"], row, color="#ff7f0e", marker="D", s=80, zorder=3,
                   label="Information pattern" if row == 0 else None)
    ax.set_yticks(range(len(results)))
    ax.set_yticklabels([r["name"].split()[0][:25] for r in results])
    ax.invert_yaxis()
    ax.set_xlabel("Position of the center of weight (bp)")
    ax.set_title("Centers of weight of each pattern")
    ax.grid(True, axis="x", linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    return fig


def write_csv(results, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sequence", "window", "start", "end", "center", "gc_percent", "information_content"])
        for r in results:
            for idx, w in enumerate(r["pattern"], 1):
                writer.writerow([r["name"].split()[0], idx, w["start"], w["end"], w["center"],
                                 round(w["gc"], 4), round(w["info"], 4)])


def main():
    parser = argparse.ArgumentParser(description="Promoter DNA patterns: C+G% and information content in sliding windows")
    parser.add_argument("fasta", nargs="?", default=None, help="FASTA file with one or more promoter sequences (default: the test sequence from the board)")
    parser.add_argument("--window", type=int, default=WINDOW, help="window length in bases (default 30)")
    parser.add_argument("--step", type=int, default=STEP, help="step between windows (default 1)")
    parser.add_argument("--csv", default=None, help="write every window to a CSV file")
    parser.add_argument("--full", action="store_true", help="print every window")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG files but do not open windows")
    args = parser.parse_args()

    if args.fasta is None:
        records = [("Board test sequence", BOARD_SEQUENCE)]
    else:
        path = args.fasta
        if not os.path.isfile(path):
            beside_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
            if os.path.isfile(beside_script):
                path = beside_script
        if not os.path.isfile(path):
            raise SystemExit(f"File not found: {args.fasta}")
        records = read_fasta(path)
        if not records:
            raise SystemExit("No FASTA records found (the file should start with '>').")

    results = []
    for name, seq in records:
        if len(seq) < args.window:
            print(f"Skipping {name[:40]}: only {len(seq)} bp, shorter than the window ({args.window}).")
            continue
        result = analyse(name, seq, args.window, args.step)
        results.append(result)
        print_summary(result, args.window)
        if len(records) == 1 or args.full:
            print_window_table(result, args.full)

    if not results:
        raise SystemExit("Nothing to analyse.")

    for i, r in enumerate(results[:MAX_PATTERN_PLOTS], 1):
        outfile = "lab8_pattern.png" if len(results) == 1 else f"lab8_pattern_{i}.png"
        draw_pattern(r, outfile)
        print(f"Saved pattern chart to {outfile}")
    draw_centers(results, "lab8_centers.png")
    print("Saved centers chart to lab8_centers.png")

    if args.csv:
        write_csv(results, args.csv)
        print(f"Saved window values to {args.csv}")
    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()