Lab 3 (variation) - Codon frequency analyzer
Reads a FASTA file, counts codons across all 6 reading frames
(3 forward + 3 on the reverse complement), prints the top N codons
with their amino acids, summarizes amino-acid usage, and plots the result.

Usage:
    python lab3_codon_frequency.py                      # uses default file below
    python lab3_codon_frequency.py my_genome.fna -n 15  # custom file / top N
    python lab3_codon_frequency.py --no-plot            # skip the matplotlib window

Requires:
    pip install matplotlib