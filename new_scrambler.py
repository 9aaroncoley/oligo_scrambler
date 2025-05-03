import random
import pandas as pd

random.seed(18)

oligo = "CAGATGCTCAACATTCCGCAGACCTCTCTGCAAGCAAAGCCCGTGGCCCCACAGGTGCCCAGCCCAGGGGGCGCCCCGGGCCAGGGTCCATACCCGTACAGCCTCTCTGAGCCAGCACCTCTCACTTTGGACACGAGCGGGAAGAATCTGACGGAGCAGAACAGCTACAGCAACATTCCTCACGAAGGGAAGCACACGCCGCTGTATGAGCGGTCCTTGCCCATCAACCCGGCCCAGAGCGGCAGCCCCAACCACGTGGATTCCGCCTACTTCCCTGGCTCTTCTACATCGTCATCTTCCGACAACGACGAGGGC"
scramble_size = 6   # we want to scramble 6 bases
shift = 2 # how many bases we want to shift over

og_seq = []
scrambled_list = [] # add to this list with stuff
id = []
bp_coordinates = []
scrambled_seq = []
hamming_collection = {}
initial_base = []
left_flanking = [] # left and right flanking sequences
right_flanking = []

id_num = -1
max_hamming = -1
for i in range(0, len(oligo) - scramble_size + 1, shift):
    id_num += 1

    replacement = oligo[i:i+scramble_size] # all the original letters
    og_seq.append(replacement) # add the chunk to oligo sequencing list
    og_temp = replacement # original sequence

    # need to calculate maximum frequency for later equation
    replacement = list(og_temp) # list of characters
    counts = {"A": replacement.count('A'), "C": replacement.count('C'), "G": replacement.count('G'), "T": replacement.count('T') }
    max_freq = max(counts.values())

    # use max letter frequency to determine which equation to use
    
    if (max_freq > (scramble_size / 2)): # if max_freq is more than half, then only half the bases can be replaced 
        max_hamming = (2 * (scramble_size - max_freq))
    else: # if max class is less than half then all bases can be replaced
        max_hamming = scramble_size
    
    at_the_maximum = False
    seen = [] # cut loops by not checking the same thing twice

    
    while (not at_the_maximum):
        random.shuffle(replacement) # swap things around
        rep_tuple = tuple(replacement)
        current_hamming = 0
        if (rep_tuple not in seen):
            seen.append(rep_tuple)
            for char1, char2, in zip(replacement, og_temp):
                if (char1 != char2):
                    current_hamming = current_hamming + 1
            if (current_hamming == max_hamming):
                at_the_maximum = True
    # we eventually get the maximum hamming

    
    
    bp_coordinates.append(str(i+1) + "-" + str(i+scramble_size)) 
    scrambled_seq.append(replacement)
    scrambled_list.append(oligo[:i] + ''.join(replacement) + oligo[i+scramble_size:])
    initial_base.append((i+1, (i+scramble_size)))
    id.append("FAM120A2P_" + str(id_num))
    left_flanking.append(oligo[:i])
    right_flanking.append(oligo[(i+scramble_size-1):])

#calculate all the AT and CG hamming distances 

hamming = [0] * len(og_seq)
at_hamming = [0] * len(og_seq)
cg_hamming = [0] * len(og_seq)
count = 0
for string1, string2 in zip(og_seq, scrambled_seq):
    for char1, char2, in zip(string1, string2):
        if (char1 != char2):
            hamming[count] = hamming[count] + 1
            
            if (char1 == 'A' or char1 == 'T'):
                at_hamming[count] = at_hamming[count] + 1
            elif (char1 == 'C' or char1 == 'G'):
                cg_hamming[count] = cg_hamming[count] + 1
    count = count + 1


output_file = pd.DataFrame()


output_file['ID'] = id
output_file['og seq'] = og_seq
output_file['scrambled seq'] = [''.join(seq) for seq in scrambled_seq]
output_file['Base Indexes'] = initial_base
output_file["Hamming Distance"] = hamming
output_file["A-T Bond Lost"] = at_hamming
output_file["C-G Bond Lost"] = cg_hamming
output_file['Left Flanking'] = left_flanking
output_file['Right Flanking'] = right_flanking
output_file['Original Oligo'] = oligo
output_file['New Oligo'] = scrambled_list
output_file_path = r"C:\Users\that9\OneDrive\Documents\Programming\Xander Scrambler\output_locations.csv"
output_file.to_csv(output_file_path, index=False)




    



