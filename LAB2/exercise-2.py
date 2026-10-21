import math

def calculate_alternative_tm(sequence, na_concentration):
    sequence = sequence.upper()
    length = len(sequence)
    
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_percent = ((g_count + c_count) / length) * 100
    
    tm = 81.5 + 16.6 * math.log10(na_concentration) + 0.41 * gc_percent - (600 / length)
    
    return tm

sequence = "TAACGCTA"
na_concentration = 0.05
tm_value = calculate_alternative_tm(sequence, na_concentration)
print(tm_value)