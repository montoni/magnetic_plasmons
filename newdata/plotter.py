import numpy as np 
import matplotlib.pyplot as plt 

lowest = np.zeros(12)
second = np.zeros(12)
for x in range(0,12):
	A = np.loadtxt('_'.join([str(x+1),'nm_eigen.txt']))
	lowest[x] = A[3]
	second[x] = A[2]

'''for x in range(10,19):
	A = np.loadtxt('_'.join([str(x+1),'nm_eigen.txt']))
	NN[x] = A[19]
	NS[x] = A[18]'''

np.savetxt('lowest.txt',lowest)
np.savetxt('second.txt',second)