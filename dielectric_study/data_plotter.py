import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

NN_nearest = np.loadtxt('NN_eps_1.txt')
NS_nearest = np.loadtxt('NS_eps_1.txt')
NN_next = np.loadtxt('../unfused_twomer/NN_corners')
NS_next = np.loadtxt('../unfused_twomer/NS_corners')
NN_next_next = np.loadtxt('../unfused_twomer/NN_onespace')
NS_next_next = np.loadtxt('../unfused_twomer/NS_onespace')

r = np.linspace(1,30,len(NN_nearest))
plt.figure(1)
plt.plot(r,NN_nearest,lw=3,label='North-North Theory')
plt.plot(r,NS_nearest,lw=3,label='North-South Theory')
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
plt.xlim([1,30])
plt.legend()
plt.savefig('NN_NS_theory.pdf')

plt.figure(2)
plt.plot(r,NN_next,lw=3,label='North-North Corner')
plt.plot(r,NS_next,lw=3,label='North-South Corner')
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
plt.xlim([1,30])
plt.legend()
plt.savefig('NN_NS_corners.pdf')

plt.figure(3)
plt.plot(r,NN_next_next,lw=3,label='North-North One Space')
plt.plot(r,NS_next_next,lw=3,label='North-South One Space')
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
plt.xlim([1,30])
plt.legend()
plt.savefig('NN_NS_spaced.pdf')
plt.show()