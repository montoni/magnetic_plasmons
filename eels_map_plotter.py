import numpy as np
import matplotlib.pyplot as plt

A = np.transpose(np.loadtxt('greybush_2017/10nm_13_EELS_map_3.52'))
AA = np.append(np.fliplr(A),A,axis=1)
eelsmap = np.append(np.flipud(AA),AA,axis=0)
plt.figure()
plt.contourf(eelsmap,cmap='hot')
plt.savefig('greybush_2017/10nm_13_aN_EELS_map.pdf')
plt.show()