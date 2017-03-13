import numpy as np
import matplotlib.pyplot as plt

for num in [1,4,7,10]:
	A = np.loadtxt('nm' + str(num) + '_twomer')
	plt.figure()
	plt.plot(A[:,0],A[:,1],A[:,0],A[:,2])
	plt.show()