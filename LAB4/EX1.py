import random
import urllib.request

def fetch_ncbi_fasta(accession="NC_001416.1"):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    try:
        with urllib.request.urlopen(url) as response:
            lines = response.read().decode('utf-8').splitlines()
            return "".join(line.strip() for line in lines if not line.startswith(">"))
    except Exception:
        return "".join(random.choices("ACGT", k=10000))

def generate_samples(sequence, num_samples=2000, min_len=100, max_len=150):
    samples = []
    seq_len = len(sequence)
    for _ in range(num_samples):
        read_len = random.randint(min_len, max_len)
        start_pos = random.randint(0, seq_len - read_len)
        samples.append(sequence[start_pos : start_pos + read_len])
    return samples

def find_max_overlap(a, b, min_overlap=15):
    start = 0
    while True:
        start = a.find(b[:min_overlap], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def greedy_assemble(reads, min_overlap=15):
    pool = list(set(reads))
    
    while len(pool) > 1:
        best_len = 0
        best_pair = None
        merged_read = None
        
        for i in range(len(pool)):
            for j in range(len(pool)):
                if i == j:
                    continue
                ov_len = find_max_overlap(pool[i], pool[j], min_overlap)
                if ov_len > best_len:
                    best_len = ov_len
                    best_pair = (i, j)
                    merged_read = pool[i] + pool[j][ov_len:]
                    
        if best_pair is None:
            break
            
        i, j = best_pair
        for index in sorted([i, j], reverse=True):
            pool.pop(index)
        pool.append(merged_read)
        
    return pool

genome = fetch_ncbi_fasta("NC_045512.2")
genome_subset = genome[:10000]

samples = generate_samples(genome_subset, num_samples=2000, min_len=100, max_len=150)
contigs = greedy_assemble(samples, min_overlap=20)

longest_contig = max(contigs, key=len)

print(f"Original sequence length: {len(genome_subset)}")
print(f"Number of samples collected: {len(samples)}")
print(f"Number of separate contigs left after assembly: {len(contigs)}")
print(f"Longest reconstructed contig length: {len(longest_contig)}")