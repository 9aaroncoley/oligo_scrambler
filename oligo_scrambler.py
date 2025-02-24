import random
import pandas as pd

random.seed(18)

# TO DO : MAXIMIZE HAMMING DISTANCE


# need to scramble some oligo and shift as necessary

# we're going tZo just create a list of the changed oligos and just keep adding to it

oligo = "CAGATGCTCAACATTCCGCAGACCTCTCTGCAAGCAAAGCCCGTGGCCCCACAGGTGCCCAGCCCAGGGGGCGCCCCGGGCCAGGGTCCATACCCGTACAGCCTCTCTGAGCCAGCACCTCTCACTTTGGACACGAGCGGGAAGAATCTGACGGAGCAGAACAGCTACAGCAACATTCCTCACGAAGGGAAGCACACGCCGCTGTATGAGCGGTCCTTGCCCATCAACCCGGCCCAGAGCGGCAGCCCCAACCACGTGGATTCCGCCTACTTCCCTGGCTCTTCTACATCGTCATCTTCCGACAACGACGAGGGC"
scrambled_list = [] # add to this list with stuff
scramble_size = 6   # we want to scramble 6 bases
shift = 2 # how many bases we want to shift over
id = []
bp_coordinates = []
og_seq = []
scrambled_seq = []
hamming_collection = {}


count = -1

for i in range(0, len(oligo) - scramble_size + 1, shift): 
    count = count + 1
    replacement = oligo[i:i+scramble_size] # has all the og letters
    og_seq.append(replacement) # add the chunk to oligo sequencing list
    og_temp = replacement # original sequence
    while True: ## need to make it so we get the max hamming distance every time
        
        
        temp = list(replacement) # list of characters
        counts = {"A": temp.count('A'), "C": temp.count('C'), "G": temp.count('G'), "T": temp.count('T') }
        non_zero_counts = [count for count in counts.values() if count != 0] # properly counts the frequency of each base
        
        if (len(set(non_zero_counts)) == 1): # only one common length, hamming max = 6
            max_hamming = 6
        else:
            max_hamming = (2 * min(non_zero_counts)) # can't find an equation, this is a good baseline
        random.shuffle(temp) # order is randomized 
        replacement = ''.join(temp)
       
        # EDGE CASE - all bases are the same
        
        if (replacement != oligo[i:i+scramble_size]): # they arent equal to each other, aka we're actually doing something
            # need max hamming
            at_diff_count = 0
            cg_diff_count = 0
            current_hamming = 0
            for string1, string2 in zip(og_temp, temp):
                for char1, char2, in zip(string1, string2):
                    if (char1 != char2):
                        current_hamming = current_hamming + 1
                        if (char1 == 'A' or char1 == 'T'):
                            at_diff_count += 1
                        elif (char1 == 'C' or char1 == 'G'):
                            cg_diff_count += 1
            if (current_hamming >= max_hamming): # if > but != 6, we wanna just try again
                if (current_hamming == 6):
                    break
                else: # >= max_hamming but not 6
                     hamming_collection.setdefault(og_temp, [])
                     if (len(hamming_collection[og_temp]) == 720): # we tried 50 times
                        # need to actually get a good oligo
                        current_max_hamming = 0
                        current_replacement = ""
                        current_max_cg_diff = 0
                        for value in hamming_collection[og_temp]:
                            if (value[1] == current_max_hamming): #same hamming, check cg
                                if (value[2] > current_max_cg_diff): # more cg differences
                                    current_max_cg_diff = value[2]
                                    current_max_hamming = value[1]
                                    current_replacement = value[0]
                            if (value[1] > current_max_hamming):
                                current_max_cg_diff = value[2]
                                current_max_hamming = value[1]
                                current_replacement = value[0]
                        replacement = current_replacement
                        break
                     temp_list = [replacement, current_hamming, at_diff_count, cg_diff_count]
                     hamming_collection.setdefault(og_temp, []).append(temp_list) # add replacement and hamming to list
                
        if (replacement[0] * scramble_size == oligo[i:i+scramble_size]): #if it's all just the same letter fr no scrambling can be done
            break
    id.append("FAM120A2P_" + str(count))
    bp_coordinates.append(str(i+1) + "-" + str(i+scramble_size)) 
    scrambled_seq.append(replacement)
    scrambled_list.append(oligo[:i] + replacement + oligo[i+scramble_size:])


output_file = pd.DataFrame()



output_file['ID'] = id
output_file['og seq'] = og_seq
output_file['scrambled seq'] = scrambled_seq

# still need to calculate the hamming distance
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

output_file["Hamming Distance"] = hamming
output_file["A-T Bond Lost"] = at_hamming
output_file["C-G Bond Lost"] = cg_hamming
output_file['Oligo'] = scrambled_list
output_file_path = "output_locations.csv"
output_file.to_csv(output_file_path, index=False)
# need id, bp scrambled, og seq segment, scrambled seq segment, and the full thing

