import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm
poynting_points = 46
mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
interaction = [[],[]]
elec = 1.60217662e-19 # regular coulombs
numPart = 13; #number of particles
me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; # does this have units?
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar*16)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
for epsb in range(1,2):
	NN = []
	NS = []
	r = 20
	for dist in np.linspace(0,4,9):
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print dist
		NN_count = 0
		NS_count = 0
		# make unit vectors in centimeters.
		rij = (dist+2)*a0
		e1 = float(1)/float(2) * rij ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * rij; #long side of 30-60-90 triangle
		Loc = [np.array([0, e1]),np.array([-e2, 2*e1]),
			   np.array([-2*e2, e1]),np.array([-2*e2, -e1]),
			   np.array([-e2, -2*e1]),np.array([0, -e1]),
			   np.array([e2, -2*e1]),np.array([2*e2, -e1]),
			   np.array([2*e2 , e1]),np.array([e2 , 2*e1]),
			   np.array([-e2,4*e1]),np.array([0,5*e1]),np.array([e2,4*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.
		center = [np.array([-e2, 0]),np.array([e2,0]),np.array([0, 3*e1])]
		'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
		H = np.zeros((2*numPart,2*numPart),dtype=complex) #initialize Hammy with zeros, twenty by twenty in this case.
		alphasp = a0**3 * (3/(epsinf+2))
		'''More constants'''

		count = 1
		#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		wsp_0 = (mie_omegas[index])*elec/hbar
		'''initialize w_0 and eigen'''
		w_mode = 0
		
		vec_NS_y = np.loadtxt('vec_mode_0')
		vec_NS_x = np.loadtxt('vec_mode_1')
		vec_NNN = np.loadtxt('vec_mode_2')
		vectors = [[vec_NS_y],[vec_NS_x],[vec_NNN]]
		coupling = 0
		count = 0
		Bfield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		Efield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		for mode in [0,2]:
			vec = np.reshape(vectors[mode],(numPart,2))
			while abs(np.real(w_mode)*hbar/elec - (np.real(wsp_0)*hbar/elec + np.real(coupling))) > 0.000001:
				if count == 0:
					w_mode = 0
					count += 1
					print count
				else:
					w_mode = np.real((wsp_0*hbar/elec + coupling) * elec/hbar)
					count += 1
					#print count
					print w_mode*hbar/elec
				wavenumber = math.sqrt(epsb)*w_mode/c
				alpha = ((alphasp**-1) - 1j*(2./3.)*wavenumber**3)**-1
				coupling = 0
				nearfield = 0
				midfield = 0
				farfield = 0
				for x in range(0,numPart):
					for y in range(x,numPart):
						if x == y:
							continue
						else:
							pass
						Rmag = math.hypot(Loc[x][0]-Loc[y][0], Loc[x][1]-Loc[y][1])
						unit_vector = (Loc[x] - Loc[y])/Rmag
						unit_dyad_term = np.dot(vec[x],vec[y])
						n_dyad_term = np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*wavenumber**2)/(Rmag)
						exponent = np.exp(1j*wavenumber*Rmag)
						coupling += -(hbar/elec)*wsp_0/epsb*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
						nearfield += -(hbar/elec) * wsp_0/epsb * r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
						midfield += -(hbar/elec) * wsp_0/epsb * -1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
						farfield += -(hbar/elec) * wsp_0/epsb * r_unit * exponent * (unit_dyad_term - n_dyad_term)
			screen_dist = 1000*a0
		
			for number in range(0,numPart):
				location = list(Loc[number])
				location.append(0)
				
				Bfield = np.empty((poynting_points,poynting_points,3),dtype=complex)
				Efield = np.empty((poynting_points,poynting_points,3),dtype=complex)
				phi_count = 0
				
				for phi in np.linspace(0,2*math.pi,poynting_points):
					theta_count = 0
					for theta in np.linspace(0,math.pi,poynting_points):
						nhat = [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]
						point = np.multiply(screen_dist,nhat)
						rmag = np.sqrt((point[0]-location[0])**2 + (point[1]-location[1])**2 + point[2]**2)
						nhat_dip = (location-point)/rmag
						Bfield[theta_count,phi_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)/rmag #
						Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)/(epsb)
						theta_count += 1
					phi_count += 1
				Bfield_total = Bfield_total + Bfield
				Efield_total = Efield_total + Efield
		poynting = np.empty((poynting_points,poynting_points),dtype=float)
		for idx1 in range(poynting_points):
			for idx2 in range(poynting_points):
				poynting[idx1,idx2] = (screen_dist**2)*np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))
		theta = np.linspace(0,math.pi,poynting_points)
		phi = np.linspace(0,2*math.pi,poynting_points)
		PHI,THETA = np.meshgrid(phi,theta)
		X = poynting * np.cos(PHI)*np.sin(THETA)
		Y = poynting * np.sin(PHI)*np.sin(THETA)
		Z = poynting * np.cos(THETA)
		norm = poynting/poynting.max()
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1, projection='3d')
		plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm),linewidth=0, antialiased=False, alpha=0.5)
		plt.show()
