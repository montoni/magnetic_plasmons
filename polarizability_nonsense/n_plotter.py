import numpy as np
import matplotlib.pyplot as plt

real_data = np.loadtxt('real_n_twomer')
imag_data = np.loadtxt('imag_n')

r = np.linspace(1,29,29)
eV = np.linspace(1,6,1001)
levels = np.linspace(-25,25,101)

plt.contourf(eV,r,real_data,levels,cmap='bwr')
plt.colorbar()

plt.show()