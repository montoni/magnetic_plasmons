import numpy as np 
import matplotlib.pyplot as plt 

A = np.loadtxt('NN.txt')
B = np.loadtxt('NS.txt')
x = np.arange(1,20)
plt.plot(x,A,'b',x,B,'purple',linewidth=3,)
#plt.axis([0.6,2.4,0,1.1])
plt.show()
