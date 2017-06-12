import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
elec = 1.60217662e-19
me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; 
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar*16)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
epsb = 1


for r in range(1,2):
	preset = {'no_coupling':1,'nearest':4,'nextnear':10}
	numPart = preset['nearest']
	radius = r*1e-7 #(in centimeters)
	rnn = 3*radius
	ex = rnn*math.sqrt(3)/2
	ey = rnn/2.
	if numPart == 1:
		loc = [0,0]
	elif numPart == 4:
		loc = [[0,0],[ex,ey],[-ex,ey],[0,-2*ey]]
	elif numPart == 10:
		loc = [[0,0],[ex,ey],[-ex,ey],[0,-2*ey],[2*ex,0],[ex,3*ey],[-ex,3*ey],[-2*ex,0],[-ex,-3*ey],[ex,3*ey]]
	else:
		pass

	'''pick a polarization direction'''
	pol = [0,0] # x-directed for now

	wsp_0 = mie_omegas[(r-1)*10]
	coupling = 0
	alpha_static = (radius**3)*(3./(epsinf+2))
	omega_system = 0
	H_int = np.zeros((numPart,numPart))
	count = 0
	while abs(np.real(omega_system)*hbar/elec - np.real(wsp_0+coupling)) > 0.00001:
		if count == 0:
			omega_system = 0
			count += 1
		else:
			omega_system = (wsp_0+coupling)*elec/hbar
			count += 1
		print count
		coupling = 0
		wavenumber = math.sqrt(epsb)*omega_system/c
		alpha_ret = ((alpha_static**-1) - 1j*(2./3.)*wavenumber**3)**-1
		for n in range(0,numPart):
			for m in range(n,numPart):
				if n == m:
					pass
				else:
					sep_vec = np.subtract(loc[n],loc[m])
					sep_mag = np.sqrt(sep_vec[0]**2 + sep_vec[1]**2)
					unit_vec = sep_vec/sep_mag
					p_dot_p = 1
					p_nn_p = np.dot(pol,unit_vec)*np.dot(unit_vec,pol)
					near = alpha_ret *sep_mag**-3
					mid = alpha_ret*wavenumber*sep_mag**-2
					far = alpha_ret*wavenumber**2*sep_mag**-1
					expon = np.exp(sep_mag*wavenumber*1j)
					coupling += -wsp_0*np.real((near * (3*p_nn_p - p_dot_p) - 1j * mid * (3*p_nn_p - p_dot_p) + far * (-p_nn_p + p_dot_p))*expon)
					H_int[n,m] = -np.real((near * (3*p_nn_p - p_dot_p) - 1j * mid * (3*p_nn_p - p_dot_p) + far * (-p_nn_p + p_dot_p))*expon)
					#print H_int[n,m]
	#w_mode = np.real(wsp_0+coupling)
	k = 1/rnn
	grid_size = 121
	omega_plus = np.zeros((grid_size,grid_size))
	omega_minus = np.zeros((grid_size,grid_size))
	i_count = 0
	for iii in np.linspace(-math.pi,math.pi,grid_size):
		i_count += 1
		j_count = 0
		for jjj in np.linspace(-math.pi,math.pi,grid_size):
			j_count += 1
			k_vector = [iii*k,jjj*k]
			loc_count = 0
			k_dot_products = 0
			for neighbor in loc:
				if neighbor == loc[0]:
					pass
				else:
					loc_count += 1
					k_dot_products += H_int[0,loc_count]*(np.cos((neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1])))	
			omega_plus[i_count-1,j_count-1] = wsp_0 + wsp_0*k_dot_products#wsp_0*np.sqrt(1 + 2*k_dot_products)
			omega_minus[i_count-1,j_count-1] = wsp_0 - wsp_0*k_dot_products#wsp_0*np.sqrt(1 - 2*k_dot_products)
	[iii,jjj] = np.meshgrid(np.linspace(-math.pi,math.pi,grid_size),np.linspace(-math.pi,math.pi,grid_size))
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(iii,jjj,omega_minus,cmap='bwr',linewidth=0)
	plt.show()


#plt.contourf(iii,jjj,omega_minus,cmap='bwr')
	#plt.colorbar()
