import numpy as np 
import matplotlib.pyplot as plt 

cross = np.loadtxt('crossing')
eps = np.arange(1,16)
plt.plot(np.sqrt(eps),cross,linewidth=3)
plt.show()


'''all_N = np.loadtxt('epsilon_1/altl_N_eigen')
alt_NS = np.loadtxt('epsilon_1/alt_NS_eigen')
out_N_in_S = np.loadtxt('epsilon_1/out_N_in_S_eigen')
dipole = np.loadtxt('epsilon_1/dipole_eigen')
nodipole = np.loadtxt('epsilon_1/nodipole_eigen')

r = np.linspace(1,30,291)
dipole = np.reshape(dipole,[291,2])
nodipole = np.reshape(nodipole,[291,2])
plt.plot(r,all_N,r,alt_NS,r,out_N_in_S,r,dipole[:,1],r,nodipole[:,1],linewidth=3)
plt.show()'''