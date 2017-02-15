import numpy as np
import matplotlib.pyplot as plt

real_data = np.loadtxt('real_n_twomer')
#imag_data = np.loadtxt('imag_n')

r = np.linspace(1,25,25)
eV = np.linspace(1,6,1001)
levels = np.linspace(-10,10,41)

plt.figure()
plt.contourf(eV,r,real_data[:][0:25],levels,cmap='bwr')
plt.colorbar()
plt.xlabel('Energy (eV)')
plt.ylabel('Scale (nm)')
plt.savefig('n_twomer.pdf')
plt.show()