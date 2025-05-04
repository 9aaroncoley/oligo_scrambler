import pandas as pd
import random

def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def bond_difference(original, scrambled):
    at_original = original.count('A') + original.count('T')
    cg_original = original.count('C') + original.count('G')
    at_scrambled = scrambled.count('A') + scrambled.count('T')
    cg_scrambled = scrambled.count('C') + scrambled.count('G')
    return abs(cg_original - cg_scrambled), abs(at_original - at_scrambled)

# Parameters
input_file = "xander_scramble_input.csv"
output_csv = "xander_scramble_output.csv"
top_n = 3  # how many scrambles per sequence

df = pd.read_csv(input_file)
data = []

for idx, row in df.iterrows():
    original = row["20nt protospacer"]
    name = row["Name"]
    
    seen = set()
    attempts = 0
    max_attempts = top_n * 10
    while len(seen) < top_n and attempts < max_attempts:
        scrambled = ''.join(c if c == 'X' else random.choice('ACGT') for c in original)
        if scrambled != original and scrambled not in seen:
            seen.add(scrambled)
        attempts += 1

    for scrambled in list(seen)[:top_n]:
        h_dist = hamming_distance(original, scrambled)
        cg_diff, at_diff = bond_difference(original, scrambled)
        data.append({
            "Original Sequence": original,
            "Scrambled Sequence": scrambled,
            "Name": name,
            "Hamming Distance": h_dist,
            "C-G Bond Lost": cg_diff,
            "A-T Bond Lost": at_diff
        })

# Create DataFrame from collected rows
output_df = pd.DataFrame(data)
output_df.to_csv(output_csv, index=False)
print(f"Output saved to {output_csv}")
