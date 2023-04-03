import numpy as np
import matplotlib.pyplot as plt
import os

fig, ax = plt.subplots()
ax.plot([1,2,4, 8,16, 32], [1,2,4, 8,16, 32], 'b-', label='Ideal')

ax.set_xlabel('Number of PEs')
ax.set_ylabel('Speedup')
ax.legend()

plt.show()
