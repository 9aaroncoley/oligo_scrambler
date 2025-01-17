import random
import pandas as pd

random.seed(18)


# need to scramble some oligo and shift as necessary

# we're going tZo just create a list of the changed oligos and just keep adding to it

test = "XXXXXXXXXXXXXXXXXXXXXXXXX"
scrambled_list = [] # add to this list with stuff
scramble_size = 6   # we want to scramble 6 bases


for i in range(0, len(test) - scramble_size + 1, 2):
    
    while True:
        
        replacement = test[i:i+scramble_size] # has all the og letters
        temp = list(replacement) # list of characters
        random.shuffle(temp) # order is randomized 
        replacement = ''.join(temp)
        
        if (replacement != test[i:i+scramble_size]): # they arent equal to each other, aka we're actually doing something
            break
        
    scrambled_list.append(test[:i] + replacement + test[i+scramble_size:])

print(scrambled_list)