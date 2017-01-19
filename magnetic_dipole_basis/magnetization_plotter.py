import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

A = np.loadtxt('chain_magnetization_eps1')
B = np.loadtxt('chain_magnetization_eps2')
C = np.loadtxt('chain_magnetization_eps3')
D = np.loadtxt('chain_magnetization_eps4')
E = np.loadtxt('chain_magnetization_eps5')
F = np.loadtxt('chain_magnetization_eps9')

r = np.linspace(1,30,291)
plt.figure(1)
plt.plot(r,A,linewidth=3,label='Epsb=1')
plt.plot(r,B,linewidth=3,label='Epsb=2')
plt.plot(r,C,linewidth=3,label='Epsb=3')
plt.plot(r,D,linewidth=3,label='Epsb=4')
plt.plot(r,E,linewidth=3,label='Epsb=5')
plt.plot(r,F,linewidth=3,label='Epsb=9')
plt.ylim([-0.1,1.1])
plt.xlim([1,30])
plt.legend(loc=2)
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Normalized Net Magnetization (unitless)')
plt.show()