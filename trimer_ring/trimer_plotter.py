import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scint

NN_theory = np.loadtxt('NN_eps_1.txt')
NS_theory = np.loadtxt('NS_eps_1.txt')
A = np.loadtxt('eigenTrimer',skiprows=1)
rad = A[:,0]
NN_BEM = A[:,1]
NS_BEM = A[:,2]

r = np.linspace(1,30/2.2,len(NN_theory))
plt.figure(1)
plt.plot(r,NN_theory,lw=3,label='All North Theory')
plt.plot(r,NS_theory[:,0],lw=3,label='North-South Theory')
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
plt.xlim([1, 15])
#plt.legend()
#plt.savefig('ring_trimer_theory.pdf')

#rad = [1,2.5,5,6,7,8,9,10]

#radnew = np.linspace(1,10,291)

#NN_smooth = spline(rad,NN_BEM,radnew)
#NS_smooth = spline(rad,NS_BEM,radnew)

#NN_smooth = scint.interp1d(rad,NN_BEM,'quadratic')
#NS_smooth = scint.interp1d(rad,NS_BEM,'quadratic')


plt.plot(rad,NN_BEM,lw=3,label='All North MNPBEM')
plt.plot(rad,NS_BEM,lw=3,label='North-South MNPBEM')
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
#plt.ylim([3.02, 3.1])
plt.legend(loc=3)
#plt.savefig('ring_trimer_MNPBEM.pdf')
plt.show()