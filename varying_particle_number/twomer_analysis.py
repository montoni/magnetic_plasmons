import math
import numpy as np
import matplotlib.pyplot as plt

NN = np.loadtxt('NN_eig.txt')
NS = np.loadtxt('NS_eig.txt')
dist = np.loadtxt('distances.txt')
maxlength = dist*1e7
c = 3e17
hbar = 1.054571726e-34
elec = 1.60217662e-19
wavelength_NN = 2*math.pi*c*hbar/(NN*elec)
wavelength_NS = 2*math.pi*c*hbar/(NS*elec)

ave = (NN+NS)/2
NN_diff = NN-ave
NS_diff = NS-ave

r = np.linspace(1,30,30)
plt.figure()
plt.plot(r,NN_diff,r,NS_diff)
plt.legend(['NN','NS'])
plt.ylabel('energy (eV)')
plt.xlabel('r_0 (nm)')
plt.savefig('scale_eig_avg.pdf')
plt.show()