import random
import pandas as pd

random.seed(18)


# need to scramble some oligo and shift as necessary

# we're going tZo just create a list of the changed oligos and just keep adding to it

test = "CAGATGCTCAACATTCCGCAGACCTCTCTGCAAGCAAAGCCCGTGGCCCCACAGGTGCCCAGCCCAGGGGGCGCCCCGGGCCAGGGTCCATACCCGTACAGCCTCTCTGAGCCAGCACCTCTCACTTTGGACACGAGCGGGAAGAATCTGACGGAGCAGAACAGCTACAGCAACATTCCTCACGAAGGGAAGCACACGCCGCTGTATGAGCGGTCCTTGCCCATCAACCCGGCCCAGAGCGGCAGCCCCAACCACGTGGATTCCGCCTACTTCCCTGGCTCTTCTACATCGTCATCTTCCGACAACGACGAGGGC"
scrambled_list = [] # add to this list with stuff
scramble_size = 6   # we want to scramble 6 bases
id = []
bp_coordinates = []
og_seq = []
scrambled_seq = []


count = -1

for i in range(0, len(test) - scramble_size + 1, 2):
    count = count + 1
    replacement = test[i:i+scramble_size] # has all the og letters
    og_seq.append(replacement)
    while True:
        
        
        temp = list(replacement) # list of characters
        random.shuffle(temp) # order is randomized 
        replacement = ''.join(temp)
        # EDGE CASE - all bases are the same

        if (replacement != test[i:i+scramble_size]): # they arent equal to each other, aka we're actually doing something
            break
        if (replacement[0] * scramble_size == test[i:i+scramble_size]): #if it's all just the same letter fr no scrambling can be done
            break
    id.append("FAM120A2P_" + str(count))
    bp_coordinates.append(str(i+1) + "-" + str(i+scramble_size)) 
    scrambled_seq.append(replacement)
    scrambled_list.append(test[:i] + replacement + test[i+scramble_size:])

#print(scrambled_seq)

output_file = pd.DataFrame()



output_file['ID'] = id
output_file['og seq'] = og_seq
output_file['scrambled seq'] = scrambled_seq

# still need to calculate the hamming distance
hamming = [0] * len(og_seq)
count = 0
for string1, string2 in zip(og_seq, scrambled_seq):
    for char1, char2, in zip(string1, string2):
        if (char1 != char2):
            hamming[count] = hamming[count] + 1
    count = count + 1

output_file["Hamming Distance"] = hamming

output_file['Oligo'] = scrambled_list
output_file_path = "output_locations.csv"
output_file.to_csv(output_file_path, index=False)
# need id, bp scrambled, og seq segment, scrambled seq segment, and the full thing

