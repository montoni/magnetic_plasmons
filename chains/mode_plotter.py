import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import pandas as pd

bem_NNN = [3.614, 3.602, 3.564, 3.494, 3.412, 3.31, 3.2]
bem_N_S = [3.612, 3.601, 3.566, 3.498, 3.404, 3.294, 3.18]
bem_NSN = [3.61, 3.6, 3.57, 3.512, 3.42, 3.292, 3.12]
bem_r = [1, 5, 10, 15, 20, 25, 30]

NNN = np.loadtxt('dielectric_NNN.txt')
NSN = np.loadtxt('dielectric_NSN.txt')
N_S = np.loadtxt('dielectric_N_S.txt')

NNN_NSN_diff = np.subtract(NNN,NSN)
NNN_N_S_diff = np.subtract(NNN,N_S)
NSN_N_S_diff = np.subtract(NSN,N_S)

NNN_int = np.loadtxt('NNN_coupling.txt')
NSN_int = np.loadtxt('NSN_coupling.txt')
N_S_int = np.loadtxt('N_S_coupling.txt')

NNN_nearmid = np.loadtxt('NNN_near_mid.txt')
NNN_far = np.loadtxt('NNN_far.txt')
NSN_nearmid = np.loadtxt('NSN_near_mid.txt')
NSN_far = np.loadtxt('NSN_far.txt')
N_S_nearmid = np.loadtxt('N_S_near_mid.txt')
N_S_far = np.loadtxt('N_S_far.txt')

r = np.linspace(.1,5,50)
epsb = np.linspace(1,3,21)

NNN_smooth = spline(bem_r,bem_NNN,r)
N_S_smooth = spline(bem_r,bem_N_S,r)
NSN_smooth = spline(bem_r,bem_NSN,r)
print 'boop'
'''plt.figure()
plt.plot(r,NNN,r,NSN,r,N_S,r,NNN_smooth,r,NSN_smooth,r,N_S_smooth,linewidth=3)
plt.legend(['NNN','NSN','N-S','BEM NNN','BEM N-S','BEM NSN'])
plt.ylabel('Energy (eV)')
plt.xlabel('Particle Radius r_0 (nm)')
#plt.show()
plt.savefig('threemer_eigenvalues_new.pdf')'''
#

plt.figure()
plt.plot(r,NNN,linewidth=3,label='NNN')
plt.plot(r,NSN,linewidth=3,label='NSN')
plt.plot(r,N_S,linewidth=3,label='N_S')
'''plt.scatter(r,NNN_nearmid, color = 'C0', marker = 'o')
plt.scatter(r,NNN_far, color = 'C0', marker = 's')
plt.scatter(r,NSN_nearmid, color = 'C1', marker = 'o')
plt.scatter(r,NSN_far, color = 'C1', marker = 's')
plt.scatter(r,N_S_nearmid, color = 'C2', marker = 'o')
plt.scatter(r,N_S_far, color = 'C2', marker = 's')'''
plt.legend()
plt.ylabel('Energy (eV)')
plt.xlabel('epsilon')
plt.show()
#plt.savefig('threemer_all_interactions.pdf')

'''plt.figure()
plt.plot(epsb,NNN_NSN_diff,epsb,NNN_N_S_diff,epsb,NSN_N_S_diff,linewidth=3)
plt.legend(['NNN-NSN','NNN-N_S','NSN-N_S'])
plt.xlabel('epsilon')
plt.ylabel('energy (eV)')
#plt.show()
plt.savefig('threemer_differences_dielectric.pdf')'''