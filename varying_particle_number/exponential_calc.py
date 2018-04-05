import math
import numpy as np
import matplotlib.pyplot as plt

omega_NN = np.loadtxt('NN_spacing_eig.txt')
omega_NS = np.loadtxt('NS_spacing_eig.txt')

c = 3e17
hbar = 1.054571726e-34
elec = 1.60217662e-19
r0 = 30
spacing = np.linspace(2*r0,12*r0,51)
wavenumber_NN = (omega_NN*elec)/(hbar*c)
wavenumber_NS = (omega_NS*elec)/(hbar*c)


NN_exp = np.exp(1j*wavenumber_NN*spacing)
NS_exp = np.exp(1j*wavenumber_NS*spacing)

ave = (NN_exp + NS_exp)/2
NN_diff = NN_exp-ave
NS_diff = NS_exp-ave

plt.figure()
plt.plot(spacing,NN_diff,spacing,NS_diff)
plt.xlabel('rnn (nm)')
plt.ylabel('diff. from average')
plt.savefig('spacing_diff_ave.pdf')
plt.show()

plt.figure()
plt.plot(spacing,NN_exp,spacing,NS_exp)
plt.xlabel('rnn (nm)')
plt.ylabel('exponential value (arb. units)')
plt.savefig('spacing_exp.pdf')
plt.show()