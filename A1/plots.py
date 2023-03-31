# Read files and plot the data
import numpy as np
import matplotlib.pyplot as plt
import os
# Read the data from the file

# Get the current working directory
dir_path = os.getcwd()

# Join the file name to the path
input_file = os.path.join(dir_path, '..\test_data\output96_1_ref.txt')
output_file = os.path.join(dir_path, 'out.txt')

input = np.loadtxt(input_file, dtype=float, delimiter=' ')
output = np.loadtxt(output_file, dtype=float, delimiter=' ')
input = input[1:]

fig, ax = plt.subplots()
ax.plot(input, 'b-', label='Reference')
ax.plot(output, 'r-', label='Result')
ax.set_xlabel('Time')
ax.set_ylabel('Amplitude')
ax.legend()
plt.show()
# plt.savefig("Fig.png")
