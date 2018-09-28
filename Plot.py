from matplotlib import pyplot as PLT
import numpy as np

with open('Output_c.txt') as f:
  v = np.loadtxt(f, delimiter=" ", dtype='float', comments="#", skiprows=1, usecols=None)
v_hist = np.ravel(v)   # 'flatten' v
fig = PLT.figure()
ax1 = fig.add_subplot(111)

n, bins, patches = ax1.hist(v_hist, bins=100, normed=1, facecolor='green')
PLT.show()