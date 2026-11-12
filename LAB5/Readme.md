# Lab 5 – Gel electrophoresis and tandem repeats

Two exercises that work on a DNA sequence of 1000–3000 nucleotides taken from NCBI.

## Files

| File | Description |
|---|---|
| `lab5_ex1_gel.py` | Exercise 1: gel electrophoresis simulation |
| `lab5_ex2_repeats.py` | Exercise 2: tandem repeat finder |
| `sequence.fasta` | Your sequence downloaded from NCBI (not included) |

## Requirements

```bash
pip install matplotlib
```

## Getting a sequence from NCBI

1. Open the NCBI Nucleotide database and find a record between 1000 and 3000 bp (for example a small plasmid, gene or mRNA).
2. Choose **Send to → Complete Record → File → FASTA** and save it as `sequence.fasta` in this folder.

If the sequence is longer than 3000 nt, the scripts use 3000 nt starting at position `--start` (default 0). If it is shorter than 1000 nt, they stop with an error.

## Exercise 1 – Gel electrophoresis

```bash
python lab5_ex1_gel.py sequence.fasta
python lab5_ex1_gel.py sequence.fasta --seed 42 --no-plot
```

1. Loads the sequence
2. Takes 10 random fragments of 100–3000 bases (maximum capped at the sequence length) and stores them in a list. Fragment lengths are all different
3. Simulates migration: short fragments travel far, long fragments stay near the wells. Distance is modeled as linear in log10(length)
4. Prints a table with start, end, length and migration distance, plus a text view of the gel
5. Saves `ex1_gel.png` showing a size ladder lane and the sample lane, with longer fragments drawn brighter

| Option | Description |
|---|---|
| `--seed N` | Random seed, to get the same fragments every run |
| `--start N` | Start position when the sequence is longer than 3000 nt |
| `--no-plot` | Save the PNG without opening a window |

## Exercise 2 – Tandem repeats

```bash
python lab5_ex2_repeats.py sequence.fasta
python lab5_ex2_repeats.py sequence.fasta --no-plot
```

Finds repeated units of 3–5 bases where the whole repeat region is 6–10 bases long.

- A repeat needs at least 2 complete, back-to-back copies of the unit (for example `GAAGAAGAA` is `GAA` x 3 = 9 bases)
- Units that are themselves repeats (`AAA`, `ATAT`) are ignored
- Each region is reported once, from its leftmost position
- Regions longer than 10 bases are counted but not listed
- Repeats with different unit lengths can overlap, so both are reported

The output is a table with start, end, unit, copies and length, followed by summary counts per unit length and the most common units. A position map is saved as `ex2_repeats.png`.

| Option | Description |
|---|---|
| `--start N` | Start position when the sequence is longer than 3000 nt |
| `--no-plot` | Save the PNG without opening a window |