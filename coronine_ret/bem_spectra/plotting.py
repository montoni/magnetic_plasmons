import matplotlib.pyplot as plt
import numpy as np

for num in range(0,40):
	spec = np.loadtxt('_'.join(['Spectrum_nm',str(num+1)]))
	eV = 1.24/spec[:,1]
	abs = spec[:,3]

	plt.plot(eV,abs,linewidth=3)
	plt.legend(['absorption'])
	plt.title(str(num))
	plt.show()
	raw_input('press enter')