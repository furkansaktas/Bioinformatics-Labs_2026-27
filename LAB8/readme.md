# Lab 8 – Promoter DNA patterns (C+G% and information content)

Slides a 30-base window along a DNA sequence, computes the C+G% and the information content of every window, plots the two curves (the "pattern"), finds the center of weight of each pattern, and plots the centers on a second chart.

## Files

| File | Description |
|---|---|
| `lab8_promoter_patterns.py` | Pattern generator |
| `promoters.fasta` | Promoter sequence(s) you want to analyze (not included) |

## Requirements

```bash
pip install matplotlib
```

## How to run

```bash
python lab8_promoter_patterns.py
python lab8_promoter_patterns.py promoters.fasta
python lab8_promoter_patterns.py promoters.fasta --csv windows.csv --no-plot
python lab8_promoter_patterns.py --window 30 --step 1 --full
```

| Option | Description |
|---|---|
| `fasta` | FASTA file with one or more promoter sequences. Without it, the test sequence from the board is used |
| `--window N` | Window length in bases (default 30) |
| `--step N` | Step between windows (default 1) |
| `--csv FILE` | Write every window (start, end, center, C+G%, information content) to a CSV file |
| `--full` | Print every window |
| `--no-plot` | Save the PNG files without opening windows |

## The steps from the board

1. **Test sequence:** `CGGACTGATCTATCT` + 27 x `A` + `CGTAGCATCTATCGATCTAT` + `CTAGCGATCTA` (73 bases). It is built in and used when no FASTA file is given
2. **Sliding window** of L = 30 bases moving 1 base at a time
3. **C+G function** (`gc_content`): percentage of G and C bases in the window
4. **Information content function** (`information_content`): for a window S, with f(a, b) = 1 when the bases are equal and 0 otherwise,

   `sigma(S) = 100 - 100 * (1 / (|S|-1)) * sum over u = 1..|S|-1 of ( (sum over i = 1..|S|-u of f(x_i, x_(i+u))) / (|S|-u) )`

   For every distance u between two bases, it takes the fraction of base pairs that are equal, averages these fractions over all distances, and subtracts the result from 100. A window of one repeated base gives 0, and a mixed window gives a higher value
5. **Plot the pattern:** `lab8_pattern.png` shows the C+G% curve and the information content curve along the sequence
6. **Center of weight:** the position where the pattern balances, `sum(position x value) / sum(value)`, using the center position of each window
7. **Centers chart:** `lab8_centers.png` shows the center of weight of the C+G pattern and of the information pattern for every sequence in the FASTA file
8. **Real promoters:** run the script on your promoter FASTA file. Use `--csv` to save the window values and compare them with the tool from the lab

## Output

- Whole-sequence values, and the window average, maximum and minimum of both measures
- The center of weight of each pattern, in bases and as a percentage of the sequence length
- A table of windows (first and last rows, or all with `--full`)
- `lab8_pattern.png` (one chart per sequence, `lab8_pattern_1.png` and so on for several, up to 4) and `lab8_centers.png`

## Notes

- Positions are 1-based. The center of window i is its start plus 15.5 for a 30-base window.
- The formula on the board has the factor 100 in an unclear place. This script treats it as a percentage, so the information content runs from 0 to 100.
- Sequences shorter than the window are skipped.