# Lab 7 – Restriction enzyme digest and gel simulation

Digests a DNA sequence (a bacterial or viral genome from NCBI) with 5 restriction enzymes, reports the cleavages and fragments for every combination of enzymes, and simulates the electrophoresis gel.

## Files

| File | Description |
|---|---|
| `lab7-restriction.py` | Digest simulator |
| `sequence.fasta` | Your genome downloaded from NCBI (not included) |

## Requirements

```bash
pip install matplotlib
```

## Getting a genome from NCBI

1. Open the NCBI Nucleotide database and search for a bacterial or viral genome (for example a bacteriophage or a virus).
2. Choose **Send to → Complete Record → File → FASTA** and save it as `sequence.fasta` in this folder.

## How to run

```bash
python lab7_restriction.py sequence.fasta
python lab7_restriction.py sequence.fasta --enzymes EcoRI HindIII
python lab7_restriction.py sequence.fasta --circular --full --no-plot
python lab7_restriction.py sequence.fasta --add-enzyme SmaI CCCGGG 3
```

| Option | Description |
|---|---|
| `fasta` | FASTA file to digest (first record is used; default `sequence.fasta`) |
| `--enzymes A B ...` | Enzymes combined in the detailed report and in the "Combined" gel lane (default: all) |
| `--add-enzyme NAME SITE CUT` | Add an enzyme: name, recognized sequence, cut position inside the site (e.g. `SmaI CCCGGG 3` cuts `CCC^GGG`) |
| `--circular` | Treat the DNA as circular (bacterial chromosome or plasmid); default is linear |
| `--full` | Print every cleavage position and every fragment |
| `--no-plot` | Save the PNG without opening a window |

## Enzymes

| Enzyme | Recognized site | Cut |
|---|---|---|
| EcoRI | GAATTC | G^AATTC |
| BamHI | GGATCC | G^GATCC |
| HindIII | AAGCTT | A^AGCTT |
| TaqI | TCGA | T^CGA |
| HaeIII | GGCC | GG^CC |

## Output

1. **Cleavages per enzyme:** a table with the recognized sequence, cut site and number of cleavages of each enzyme
2. **All combinations:** every combination of the enzymes (31 for 5 enzymes) with the number of cuts, number of fragments, and smallest and largest fragment
3. **Detailed digest** of the selected combination: cleavage positions per enzyme, then a fragment table with start, end and length. The fragment lengths add up to the sequence length
4. **Gel image** `lab7_gel.png`: a size ladder, one lane per enzyme, and a "Combined" lane for the selected combination

## How it works

- Cleavage positions are the number of bases before the cut, so a position of 1,161 means the DNA is cut after base 1,161.
- A linear molecule with n cuts gives n + 1 fragments. A circular one gives n fragments, and the last fragment wraps around the end of the sequence.
- When two enzymes cut at the same position, the cut is counted once in the combination.
- The five enzymes are palindromic, so searching one strand finds every site. Added enzymes with non-palindromic sites are searched on both strands.
- The gel uses the same model as Lab 5: distance travelled is linear in log10(fragment length). Equal-sized fragments share one band, and band brightness grows with the amount of DNA in it. Fragments under 10 bp are drawn at the bottom edge.