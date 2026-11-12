import argparse
import random
from dataclasses import dataclass

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SEQ_MIN, SEQ_MAX = 200, 400
HOST_MIN, HOST_MAX = 140, 260
N_ELEMENTS_MIN, N_ELEMENTS_MAX = 3, 4

IR_SIM_MIN, IR_SIM_MAX = 8, 10
BODY_MIN, BODY_MAX = 15, 40
TSD_LEN = 5

IR_DETECT_MIN = 7

PAIR = dict(zip("ACGT", "TGCA"))
COLOURS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]


@dataclass
class Element:
    tsd_start: int
    start: int
    end: int
    tsd_end: int
    ir_len: int
    tsd: str


def revcomp(seq):
    return "".join(PAIR[b] for b in reversed(seq))


def random_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def protected_positions(elements, length):
    mask = [False] * length
    for e in elements:
        for lo, hi in ((e.tsd_start, e.start + e.ir_len), (e.end - e.ir_len, e.tsd_end)):
            for x in range(lo, hi):
                mask[x] = True
    return mask


def candidate_sites(length, elements, nested):
    blocked = protected_positions(elements, length)
    sites = []
    for p in range(length - TSD_LEN + 1):
        if any(blocked[p:p + TSD_LEN]):
            continue
        inside_body = any(e.start + e.ir_len <= p and p + TSD_LEN <= e.end - e.ir_len for e in elements)
        touches_element = any(p < e.end and p + TSD_LEN > e.start for e in elements)
        if nested and inside_body:
            sites.append(p)
        elif not nested and not touches_element:
            sites.append(p)
    return sites


def insert_element(seq, elements, rng, nested):
    sites = candidate_sites(len(seq), elements, nested)
    if not sites:
        return None
    p = rng.choice(sites)
    q = p + TSD_LEN
    ir_len = rng.randint(IR_SIM_MIN, IR_SIM_MAX)
    left_ir = random_dna(ir_len, rng)
    core = left_ir + random_dna(rng.randint(BODY_MIN, BODY_MAX), rng) + revcomp(left_ir)
    shift = len(core) + TSD_LEN
    new_seq = seq[:q] + core + seq[p:q] + seq[q:]
    for e in elements:
        if e.tsd_start >= q:
            e.tsd_start += shift
        if e.start >= q:
            e.start += shift
        if e.end > q:
            e.end += shift
        if e.tsd_end > q:
            e.tsd_end += shift
    elements.append(Element(p, q, q + len(core), q + len(core) + TSD_LEN, ir_len, seq[p:q]))
    return new_seq


def simulate(rng, n_elements, n_nested):
    while True:
        seq = random_dna(rng.randint(HOST_MIN, HOST_MAX), rng)
        nested_slots = set(rng.sample(range(1, n_elements), min(n_nested, n_elements - 1)))
        elements = []
        for k in range(n_elements):
            seq = insert_element(seq, elements, rng, k in nested_slots)
            if seq is None:
                break
        if seq is not None and SEQ_MIN <= len(seq) <= SEQ_MAX:
            return seq, sorted(elements, key=lambda e: e.start)


def find_transposons(seq):
    n = len(seq)
    found = []
    for i in range(TSD_LEN, n - 2 * IR_DETECT_MIN - TSD_LEN + 1):
        left_tsd = seq[i - TSD_LEN:i]
        for j in range(i + 2 * IR_DETECT_MIN, n - TSD_LEN + 1):
            if seq[j:j + TSD_LEN] != left_tsd:
                continue
            ir_len, limit = 0, (j - i) // 2
            while ir_len < limit and PAIR[seq[j - 1 - ir_len]] == seq[i + ir_len]:
                ir_len += 1
            if ir_len >= IR_DETECT_MIN:
                found.append({"start": i + 1, "end": j, "ir_len": ir_len, "tsd": left_tsd})
    return sorted(found, key=lambda f: (f["start"], f["end"]))


def nesting_depth(spans):
    return [sum(1 for o in spans if o is not s and o[0] < s[0] and s[1] < o[1]) for s in spans]


def container_of(span, spans):
    parents = [i for i, o in enumerate(spans) if o[0] < span[0] and span[1] < o[1]]
    return min(parents, key=lambda i: spans[i][1] - spans[i][0]) if parents else None


def print_annotated(seq, elements):
    marks = ["."] * len(seq)
    for e in sorted(elements, key=lambda e: e.start - e.end):
        label = str(elements.index(e) + 1)
        for x in range(e.start, e.end):
            marks[x] = label
    for e in elements:
        for x in list(range(e.tsd_start, e.start)) + list(range(e.end, e.tsd_end)):
            marks[x] = "="
    print("\nSequence (1-60 per line). Marks: . host, 1-4 transposon number, = target-site duplication")
    for i in range(0, len(seq), 60):
        print(f"{i + 1:>4}  {seq[i:i + 60]}")
        print(f"      {''.join(marks[i:i + 60])}")


def print_true_table(elements):
    spans = [(e.start + 1, e.end) for e in elements]
    print(f"\n{'ID':<5}{'Start':>7}{'End':>7}{'Length':>8}{'IR':>5}  {'TSD':<7}Inside")
    print("-" * 50)
    for k, (e, span) in enumerate(zip(elements, spans), 1):
        parent = container_of(span, spans)
        inside = f"T{parent + 1}" if parent is not None else "-"
        print(f"T{k:<4}{span[0]:>7}{span[1]:>7}{span[1] - span[0] + 1:>8}{e.ir_len:>5}  "
              f"{e.tsd:<7}{inside}")


def print_detected_table(found, truth_spans):
    spans = [(f["start"], f["end"]) for f in found]
    print(f"\n{'ID':<5}{'Start':>7}{'End':>7}{'Length':>8}{'IR':>5}  {'TSD':<7}{'Inside':<8}Result")
    print("-" * 62)
    for k, (f, span) in enumerate(zip(found, spans), 1):
        parent = container_of(span, spans)
        inside = f"D{parent + 1}" if parent is not None else "-"
        result = f"matches T{truth_spans.index(span) + 1}" if span in truth_spans else "not in simulation"
        print(f"D{k:<4}{f['start']:>7}{f['end']:>7}{f['end'] - f['start'] + 1:>8}{f['ir_len']:>5}  "
              f"{f['tsd']:<7}{inside:<8}{result}")


def draw_tracks(seq, truth_spans, found, outfile):
    found_spans = [(f["start"], f["end"]) for f in found]
    fig, axes = plt.subplots(2, 1, figsize=(11, 5), sharex=True)
    for ax, spans, title in ((axes[0], truth_spans, "Simulated transposons"),
                             (axes[1], found_spans, "Detected transposons")):
        depths = nesting_depth(spans)
        ax.add_patch(Rectangle((0.5, -0.45), len(seq), 0.25, color="#cccccc"))
        for k, (span, depth) in enumerate(zip(spans, depths)):
            matched = span in truth_spans
            colour = COLOURS[truth_spans.index(span) % len(COLOURS)] if matched else "#d62728"
            ax.add_patch(Rectangle((span[0] - 0.5, depth * 0.8), span[1] - span[0] + 1, 0.6,
                                   color=colour, alpha=0.85))
            name = f"T{truth_spans.index(span) + 1}" if matched else f"D{k + 1}?"
            ax.text((span[0] + span[1]) / 2, depth * 0.8 + 0.3, f"{name}\n{span[0]}-{span[1]}",
                    ha="center", va="center", color="white", fontsize=8)
        ax.set_ylim(-0.6, 0.8 * (max(depths) if depths else 0) + 1.0)
        ax.set_yticks([])
        ax.set_title(title, loc="left", fontsize=10)
    axes[1].set_xlim(0, len(seq) + 1)
    axes[1].set_xlabel("Position in sequence (nt)")
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Transposable element simulation and detection")
    parser.add_argument("--seed", type=int, default=None, help="random seed for a reproducible sequence")
    parser.add_argument("--elements", type=int, choices=[3, 4], default=None, help="number of transposons (default: random 3-4)")
    parser.add_argument("--nested", type=int, choices=[0, 1, 2, 3], default=1, help="how many transposons are inserted inside another one")
    parser.add_argument("--no-plot", action="store_true", help="save the PNG but do not open a window")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    n_elements = args.elements or rng.randint(N_ELEMENTS_MIN, N_ELEMENTS_MAX)
    seq, elements = simulate(rng, n_elements, args.nested)

    truth_spans = [(e.start + 1, e.end) for e in elements]
    print(f"Artificial DNA: {len(seq)} nt with {len(elements)} transposable elements "
          f"({sum(d > 0 for d in nesting_depth(truth_spans))} inserted inside another element)")
    print_annotated(seq, elements)
    print("\nSimulated (true) positions:")
    print_true_table(elements)

    found = find_transposons(seq)
    print("\nDetected positions:")
    print_detected_table(found, truth_spans)

    found_spans = [(f["start"], f["end"]) for f in found]
    hits = sum(span in found_spans for span in truth_spans)
    extra = sum(span not in truth_spans for span in found_spans)
    print(f"\nFound {hits} of {len(truth_spans)} transposons, {extra} false positive(s)")

    draw_tracks(seq, truth_spans, found, "lab6_transposons.png")
    print("Saved plot to lab6_transposons.png")
    if not args.no_plot:
        plt.show()


if __name__ == "__main__":
    main()