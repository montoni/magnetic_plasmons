import numpy as np
import matplotlib.pyplot as plt

for num in [10]:
	A = np.loadtxt('new_spacing_bem_data/nm' + str(num) + '_twomer')
	plt.figure()
	plt.plot(A[:,0],A[:,1],A[:,0],A[:,2])
	plt.legend(['North-South','North-North'])
	plt.show()