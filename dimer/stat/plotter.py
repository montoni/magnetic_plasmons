import numpy as np 
import matplotlib.pyplot as plt 

lowest = np.zeros(41)
second = np.zeros(41)
for x in range(0,41):
	A = np.loadtxt('_'.join([str(float(x+1)/4),'nm_eigen.txt']))
	lowest[x] = A[3]
	second[x] = A[2]

np.savetxt('lowest.txt',lowest)
np.savetxt('second.txt',second)