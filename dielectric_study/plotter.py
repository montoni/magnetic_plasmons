import numpy as np 
import matplotlib.pyplot as plt 

one = np.zeros(40)
two = np.zeros(40)
for epsb in range (1,17):
	for x in range(0,40):
		A = np.loadtxt('/'.join(['_'.join(['diel',str(epsb)]),'_'.join(['epsb',str(epsb),str(x+1),'nm_eigen.txt'])]))
		if x == 0:
			one[x] = A[19]
			two[x] = A[18]
		elif x == 1:
			one[x] = A[19]
			two[x] = A[18]
		else:
			if two[x-1]-A[18] < two[x-2]-two[x-1]:
				one[x] = A[18]
				two[x] = A[19]
			else: 
				one[x] = A[19]
				two[x] = A[18]

		
	r = np.arange(1,41)
	plt.plot(r,one,r,two,linewidth=3)
	plt.show()
