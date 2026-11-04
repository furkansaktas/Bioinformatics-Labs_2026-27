# Lab 1 (variation) - Nucleotide and k-mer frequency analysis Counts single nucleotides, GC content, and 2-mer / 3-mer frequencies
# for a short DNA sequence, then draws simple text bar charts.

from collections import Counter

dna = "GGATCCATGAAGCTTACGTAGCTA"
print("DNA Sequence:", dna)
print("Length:", len(dna))

print("\nSingle Nucleotide Frequency:")
base_counts = Counter(dna)
for base in "ATCG":
    count = base_counts.get(base, 0)
    print(f"{base}: {count} ({count / len(dna) * 100:.1f}%)")

gc = (base_counts["G"] + base_counts["C"]) / len(dna) * 100
print(f"\nGC content: {gc:.1f}%")

def kmer_frequency(sequence, k):
    kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]
    counts = Counter(kmers)
    # sort by count (high -> low), then alphabetically for ties
    return dict(sorted(counts.items(), key=lambda x: (-x[1], x[0])))


def show_chart(title, data):
    print(f"\n{title}")
    for kmer, count in data.items():
        bar = "#" * count
        print(f"{kmer}: {bar} ({count})")


freq2 = kmer_frequency(dna, 2)
freq3 = kmer_frequency(dna, 3)

show_chart("2-Letter Combination Frequency:", freq2)
show_chart("3-Letter Combination Frequency:", freq3)

top2 = next(iter(freq2.items()))
top3 = next(iter(freq3.items()))
print(f"\nMost frequent 2-mer: {top2[0]} ({top2[1]}x)")
print(f"Most frequent 3-mer: {top3[0]} ({top3[1]}x)")