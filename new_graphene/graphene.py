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
stupid_angle = np.arcsin(math.sqrt(1./3.))

for r in range(8,9):
	preset = {'no_coupling':1,'nearest':4,'nextnear':10,'nextnext':13}
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
		loc = [[0,0],[ex,ey],[-ex,ey],[0,-2*ey],[2*ex,0],[ex,3*ey],[-ex,3*ey],[-2*ex,0],[-ex,-3*ey],[ex,-3*ey]]
	elif numPart == 13:
		loc = [[0,0],[ex,ey],[-ex,ey],[0,-2*ey],[2*ex,0],[ex,3*ey],[-ex,3*ey],[-2*ex,0],[-ex,-3*ey],[ex,-3*ey],[2*ex,ey],[-2*ex,ey],[0,-4*ey]]
	else:
		pass

	'''pick a polarization direction'''
	pol = [0,math.sqrt(1./3.)] # z-directed for now

	wsp_0 = mie_omegas[(r-1)*10]
	coupling = 0
	alpha_static = (radius**3)*(3./(epsinf+2))
	omega_system = 0
	H_int = np.zeros((numPart,numPart))
	near = np.zeros((numPart,numPart))
	mid = np.zeros((numPart,numPart))
	far = np.zeros((numPart,numPart))
	p_dot_p = np.zeros((numPart,numPart))
	p_nn_p = np.zeros((numPart,numPart))
	count = 0
	while abs(np.real(omega_system)*hbar/elec - np.real(wsp_0+coupling)) > 0.000001:
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
					p_dot_p[n,m] = 1#np.dot(pol,pol)
					p_nn_p[n,m] = np.dot(pol,unit_vec)*np.dot(unit_vec,pol)
					near[n,m] = np.real(alpha_ret *sep_mag**-3)
					mid[n,m] = np.real(alpha_ret*wavenumber*sep_mag**-2)
					far[n,m] = np.real(alpha_ret*wavenumber**2*sep_mag**-1)
					expon = np.exp(sep_mag*wavenumber*1j)
					coupling += -wsp_0*np.real((near[n,m] * (3*p_nn_p[n,m] - p_dot_p[n,m]) - 1j * mid[n,m] * (3*p_nn_p[n,m] - p_dot_p[n,m]) + far[n,m] * (-p_nn_p[n,m] + p_dot_p[n,m]))*expon)
					H_int[n,m] = -np.real((near[n,m] * (3*p_nn_p[n,m] - p_dot_p[n,m]) - 1j * mid[n,m] * (3*p_nn_p[n,m] - p_dot_p[n,m]) + far[n,m] * (-p_nn_p[n,m] + p_dot_p[n,m]))*expon)
					#print H_int[n,m]
	#w_mode = np.real(wsp_0+coupling)
	k = 1/rnn
	grid_size = 251
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
			k_dot_products_near = 0
			k_dot_products_mid = 0
			k_dot_products_far = 0
			near_tot = 0
			mid_tot = 0
			far_tot = 0
			for neighbor in loc:
				if neighbor == loc[0]:
					pass
				else:
					loc_count += 1
<<<<<<< HEAD
					k_dot_products_near += (3*p_nn_p[0,loc_count] - p_dot_p[0,loc_count])*(np.exp(1j*(neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1])))
					k_dot_products_mid += (3*p_nn_p[0,loc_count] - p_dot_p[0,loc_count])*(np.exp(1j*(neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1])))
					k_dot_products_far += (-p_nn_p[0,loc_count] + p_dot_p[0,loc_count])*(np.exp(1j*(neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1])))
					near_tot += near[0,loc_count]
					mid_tot += mid[0,loc_count]
					far_tot += far[0,loc_count]
			k_dot_products = near_tot * np.absolute(k_dot_products_near) + mid_tot * np.absolute(k_dot_products_mid) + far_tot * np.absolute(k_dot_products_far)
			omega_plus[i_count-1,j_count-1] = np.sqrt(1 + 1*(k_dot_products))#wsp_0*np.sqrt(1 + 2*k_dot_products)
			omega_minus[i_count-1,j_count-1] = np.sqrt(1 - 1*(k_dot_products))#wsp_0*np.sqrt(1 - 2*k_dot_products)
=======
<<<<<<< HEAD
					k_dot_products += H_int[0,loc_count]*(np.cos((neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1])))	
=======
					k_dot_products += H_int[0,loc_count]*np.sqrt(np.square(np.cos((neighbor[0]*k_vector[0] + neighbor[1]*k_vector[1]))))
>>>>>>> 23f1ea348cb13ea4d3c878b1d00b348a5a2b4147
			omega_plus[i_count-1,j_count-1] = wsp_0 + wsp_0*k_dot_products#wsp_0*np.sqrt(1 + 2*k_dot_products)
			omega_minus[i_count-1,j_count-1] = wsp_0 - wsp_0*k_dot_products#wsp_0*np.sqrt(1 - 2*k_dot_products)
>>>>>>> e94e67c4331899ceaee706b1b970aa0f0b26448f
	[iii,jjj] = np.meshgrid(np.linspace(-math.pi,math.pi,grid_size),np.linspace(-math.pi,math.pi,grid_size))
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
<<<<<<< HEAD
	ax.plot_surface(iii,jjj,omega_minus,cmap='bwr',linewidth=0)
=======
	ax.plot_surface(iii,jjj,omega_minus,color='blue')
	ax.plot_surface(iii,jjj,omega_plus,color='red')
>>>>>>> 23f1ea348cb13ea4d3c878b1d00b348a5a2b4147
	plt.show()


#plt.contourf(iii,jjj,omega_minus,cmap='bwr')
	#plt.colorbar()
