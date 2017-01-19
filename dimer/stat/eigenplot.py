import numpy as np 
import matplotlib.pyplot as plt 

A = np.loadtxt('lowest.txt')
B = np.loadtxt('second.txt')
x = np.arange(1,42)
plt.plot(x,A,'b',x,B,'purple',linewidth=3,)
plt.xlim([0.0,42])
plt.show()
