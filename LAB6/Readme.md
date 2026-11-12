# Lab 6 – Transposable elements: simulation and detection

Builds an artificial DNA sequence of 200–400 bases containing 3–4 transposable elements, then detects the start and end position of each element using only the final sequence.

## Files

| File | Description |
|---|---|
| `lab6_transposons.py` | Simulator and detector |

## Requirements

```bash
pip install matplotlib
```

## How to run

```bash
python lab6_transposons.py
python lab6_transposons.py --seed 21 --elements 4 --nested 2
python lab6_transposons.py --no-plot
```

| Option | Description |
|---|---|
| `--seed N` | Random seed, to get the same sequence every run |
| `--elements 3\|4` | Number of transposons (default: random 3 or 4) |
| `--nested N` | How many transposons are inserted inside another transposon, 0–3 (default: 1) |
| `--no-plot` | Save the PNG without opening a window |

## How a transposon is modeled

Each element follows the scheme from the lab board:

1. The transposon is `IR + body + reverse complement of IR`. The two ends are **inverted repeats** (8–10 bases) and the body is 15–40 random bases.
2. It inserts into a 5-base **target site** of the host DNA. The site is copied, so the same 5 bases appear on both sides of the element. These are the **direct repeats** (target-site duplication).

```
host:     ...TAC ATGCA CAG...
result:   ...TAC ATGCA [IR--body--IR'] ATGCA CAG...
```

A transposon can insert into the body of an earlier one, which creates the intersecting cases (a transposon inside a transposon). Insertions never fall inside the inverted repeats or target-site duplications of an existing element, so every element stays intact. The simulator builds the sequence step by step and updates all earlier positions after each insertion, so the true positions are known exactly.

## How the detector works

The detector sees only the DNA string. For every pair of positions (i, j) it checks:

1. The 5 bases just before `i` equal the 5 bases just after `j` (direct repeats).
2. The bases starting at `i` are the reverse complement of the bases ending at `j` (inverted repeats), for at least 7 bases.

A pair that passes both checks is reported as a transposon, from the first base of the left inverted repeat (`start`) to the last base of the right inverted repeat (`end`). Positions are 1-based and inclusive. Nested elements are found independently because each one has its own target-site duplication.

## Output

- The sequence with a marker line (`.` host, `1`–`4` transposon number, `=` target-site duplication)
- A table of the true positions (start, end, length, IR length, TSD, which element it is inside)
- A table of the detected positions, with a note saying which true element each one matches
- A summary line, for example `Found 4 of 4 transposons, 0 false positive(s)`
- `lab6_transposons.png`: simulated and detected elements on two tracks, with nested elements drawn above their container

## Notes

- In a test of 6,400 random sequences with 0–3 nested insertions, every inserted element was found. About 0.15% of the sequences had one extra match caused by chance (random DNA can contain a matching inverted repeat and target-site pair). Such a hit is labeled `not in simulation` in the table.
- Insertions inside the inverted repeats of another element are not simulated, because they break the IR match that the detector relies on.
- All settings (sequence length, IR length, body length, target-site length) are constants at the top of the script.