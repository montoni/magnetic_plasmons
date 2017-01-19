import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scint

'''A = np.loadtxt('nm1_three_chain')
B = np.loadtxt('nm4_chain_three')
C = np.loadtxt('nm8_chain_three')
D = np.loadtxt('nm10_chain_three')'''

A = np.loadtxt('three_chain_eigen',skiprows=1)
NNN = np.loadtxt('../NNN.txt')
N_S = np.loadtxt('../N_S.txt')
NSN = np.loadtxt('../NSN.txt')

rad = A[:,0]
NSN_BEM = A[:,1]
N_S_BEM = A[:,2]
NNN_BEM = A[:,3]
radnew = np.linspace(1,10,291)

NNN_smooth = scint.interp1d(rad,NNN_BEM,'quadratic')
N_S_smooth = scint.interp1d(rad,N_S_BEM,'quadratic')
NSN_smooth = scint.interp1d(rad,NSN_BEM,'quadratic')

r = np.linspace(1,30,291)
plt.figure(1)
plt.plot(r,NNN,lw=3,label='NNN theory')
plt.plot(r,N_S,lw=3,label='N_S theory')
plt.plot(r,NSN,lw=3,label='NSN theory')
plt.xlim(1,30)
plt.xlabel('Radius (nm)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.savefig('three_chain_theory.pdf')
plt.figure(2)
plt.plot(radnew,NNN_smooth(radnew),lw=3,label='NNN bem')
plt.plot(radnew,N_S_smooth(radnew),lw=3,label='N_S bem')
plt.plot(radnew,NSN_smooth(radnew),lw=3,label='NSN bem')
plt.ylim([3.02,3.10])
plt.xlabel('Radius (nm)')
plt.ylabel('Energy (eV)')
plt.legend()
plt.savefig('three_chain_mnpbem.pdf')
plt.show()

'''plt.figure(1)
plt.plot(A[:,0],A[:,1],label = 'left')
plt.plot(A[:,0],A[:,2],label = 'top')
plt.legend()
plt.figure(2)
plt.plot(B[:,0],B[:,1],label = 'left')
plt.plot(B[:,0],B[:,2],label = 'top')
plt.legend()
plt.figure(3)
plt.plot(C[:,0],C[:,1],label = 'left')
plt.plot(C[:,0],C[:,2],label = 'top')
plt.legend()
plt.figure(4)
plt.plot(D[:,0],D[:,1],label = 'left')
plt.plot(D[:,0],D[:,2],label = 'top')
plt.legend()
plt.show()'''