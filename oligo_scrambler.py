import random
import pandas as pd

random.seed(18)


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


count = -1

for i in range(0, len(oligo) - scramble_size + 1, shift): 
    count = count + 1
    replacement = oligo[i:i+scramble_size] # has all the og letters
    og_seq.append(replacement)
    while True:
        
        
        temp = list(replacement) # list of characters
        random.shuffle(temp) # order is randomized 
        replacement = ''.join(temp)
       
        # EDGE CASE - all bases are the same

        if (replacement != oligo[i:i+scramble_size]): # they arent equal to each other, aka we're actually doing something
            break
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

