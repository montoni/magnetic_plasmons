import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

#index = []
collinear = []
anticollinear = []
parallel = []
antiparallel = []
count = 0
for r in range(1,41):
	bound = 0
	spectrum = np.loadtxt(''.join(['radius_',''.join([str(r),'nm'])]),skiprows=1)
	eV = spectrum[:,0]
	#print spectrum[0,1]
	#print eV
	'''plt.figure
	plt.plot(eV,spectrum[:,2],linewidth=3)
	plt.show()'''
	for i in enumerate(spectrum[:,2]):
		if i[0] < len(spectrum[:,2]):
			if spectrum[i[0],2] > bound:
				bound = spectrum[i[0],2]
				if spectrum[i[0]+1,2] < spectrum[i[0],2]:
					count = count+1
					bound = spectrum[i[0]+1,2]
		if count == 2:
			index = i[0]
			break	
	antiparallel.append(eV[index])
#,'top','middle'
r = np.arange(1,41)
plt.figure
plt.plot(r,antiparallel,linewidth=3)
plt.show()
np.savetxt('anti_para_eig',antiparallel)