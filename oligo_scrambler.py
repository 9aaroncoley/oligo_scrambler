import random


# need to scramble some oligo and shift as necessary

# we're going tZo just create a list of the changed oligos and just keep adding to it

test = "CAGATGCTCAACATTCCGCAGACCTCTCTGCAAGCAAAGCCCGTGGCCCCACAGGTGCCCAGCCCAGGGGGCGCCCCGGGCCAGGGTCCATACCCGTACAGCCTCTCTGAGCCAGCACCTCTCACTTTGGACACGAGCGGGAAGAATCTGACGGAGCAGAACAGCTACAGCAACATTCCTCACGAAGGGAAGCACACGCCGCTGTATGAGCGGTCCTTGCCCATCAACCCGGCCCAGAGCGGCAGCCCCAACCACGTGGATTCCGCCTACTTCCCTGGCTCTTCTACATCGTCATCTTCCGACAACGACGAGGGC"
scrambled_list = [] # add to this list with stuff
scramble_size = 6   # we want to scramble 6 bases
num_to_base = str.maketrans("1234", "ACGT")

for i in range(0, len(test) - scramble_size + 1, 2):
    
    while True:
        
        random_number = str(int(''.join(str(random.randint(1,4)) for num in range(scramble_size))))
        replacement = random_number.translate(num_to_base)
        
        if (replacement != test[i:i+scramble_size]): # they arent equal to each other, aka we're actually doing something
            break
        
    scrambled_list.append(test[:i] + replacement + test[i+scramble_size:])

print(scrambled_list[0])