import random
# specify the filename and the number of integers to generate
num_integers = 1000000
filename = f"{str(num_integers)}.txt"

# generate random integers and write them to the file
with open(filename, 'w') as file:
    file.write(str(num_integers) + " ")
    for i in range(num_integers):
        # generate a random integer between 1 and 100
        rand_int = random.randint(1, num_integers/10)
        # write the integer to the file followed by a space
        file.write(str(rand_int) + " ")
