import numpy as np
import matplotlib.pyplot as plt

strong = [0.706233, 0.379075, 0.602186, 0.694319, 0.086181, 0.091930]
strong = [1.0/x for x in strong]
weak = [0.708979, 0.724634, 1.129856, 1.477219, 1.282406, 2.381665]
weak = [x/weak[0] for x in weak]
weaktest = [0.007464, 0.007862, 0.024981, 0.061325]
weaktest = [x/weaktest[0] for x in weaktest]
fig1, ax1 = plt.subplots()
ax1.plot([1,2,4,8,16,32], [1,2,4, 8,16, 32], 'b--', label='Ideal')
ax1.plot([1,2,4,8,16,32], strong, 'r-', label='Strong Scaling')
ax1.set_xlabel('Number of PEs')
ax1.set_ylabel('Speedup')
ax1.set_xticks([1,2,4, 8,16, 32])
ax1.legend()
ax1.set_title("Strong Scaling")
fig1.savefig("strong.png")

fig, ax = plt.subplots()
ax.plot([1,2,4, 8,16, 32], [1,1,1, 1,1, 1], 'b--', label='Ideal')
ax.plot([1,2,4, 8,16, 32], weak, 'r-', label='Weak Scaling')
ax.set_xlabel('Number of PEs')
ax.set_ylabel('Scaled time')
ax.set_xticks([1,2,4, 8,16, 32])
ax.legend()
ax.set_title("Weak Scaling")
fig.savefig("weak.png")

fig3, ax3 = plt.subplots()
ax3.plot([1,2,4, 8], [1,1,1, 1], 'b--', label='Ideal')
ax3.plot([1,2,4, 8], weaktest, 'r-', label='Weak Scaling')
ax3.set_xlabel('Number of PEs')
ax3.set_ylabel('Scaled time')
ax3.set_xticks([1,2,4, 8])
ax3.legend()
ax3.set_title("Weak Scaling")
fig3.savefig("weaktest.png")
plt.show()
