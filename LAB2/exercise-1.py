import matplotlib.pyplot as plt

def plot_frequency(sequence):
    sequence = sequence.upper()
    bases = ['A', 'C', 'G', 'T']
    counts = [sequence.count(base) for base in bases]
    
    plt.bar(bases, counts, color=['blue', 'orange', 'green', 'red'])
    plt.title(f'Nucleotide Frequency in S = {sequence}')
    plt.xlabel('Base')
    plt.ylabel('Frequency')
    plt.show()

s = "TAACGCTA"
plot_frequency(s)