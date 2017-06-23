import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
nearest = np.loadtxt('../varying_particle_number/nearest_magnetic_coupling.txt')
#ns_near = np.loadtxt('../dielectric_study/NS_right.txt')
#nn_sec = np.loadtxt('../unfused_twomer/NN_sec_near')
#ns_sec = np.loadtxt('../unfused_twomer/NS_sec_near')
#nn_third = np.loadtxt('../unfused_twomer/NN_third_near')
#ns_third = np.loadtxt('../unfused_twomer/NS_third_near')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
all_N = []
alt_NS = []
dipole =[]
nodipole =[]
excited = []

for r in [1,15,30]:
	elec = 1.60217662e-19 # regular coulombs
	numPart = 6 #number of particles
	a0 = .1*r*10**-7 #sphere radius in cm
	mie_index = (r*10)-10
	index = (r*10)-10

	''' now determine geometry.'''
	''' begin with electric dipole on the six-member ring '''
	# make unit vectors in centimeters.
	rij = 3*a0
	e1 = float(1)/float(2) * (3*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (3*a0); #long side of 30-60-90 triangle
	m_sep = 2*e2
	Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),
		   np.array([0,-2*e1]),np.array([e2,-e1]),np.array([e2,e1])]

	Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

	Q_mag_2 = [0,0,-1]

	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = np.zeros((2*numPart,2*numPart)) #initialize Hammy with zeros, twenty by twenty in this case.

	'''More constants'''

	me = 9.10938291e-28; # electron mass in g
	ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
	hbar = 1.054571726e-34; # modified Planck in J*s
	c = 2.99e10; # speed of light in cm/s
	eps0 = 8.85418782e-12; # permittivity of free space
	epsb = 1; # background dielectric - right now we're in vacuum
	epsinf = 3.77; # does this have units?
	'''Properties for silver.'''
	Eplasma = 1.46599161e-18; # J
	gamma = 0.05*elec/(hbar*16)
	wplasma = Eplasma/hbar; # plasma frequency (rad/s)
	wsp_0 = (mie_omegas[mie_index])*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = wsp_0
	eigen = np.zeros(2*numPart)
	count = 1

	while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-1])) > 0.00000001:
		if count == 1:
			w_0 = 0*elec/hbar
			count = count + 1
			wsp = wsp_0
		else:
			count = count + 1
			w_0 = eigen[2*numPart-1]*elec/hbar
			#if w_0 > wsp:
			#	w_0 = 3.5*elec/hbar
		print eigen[2*numPart-1]
		wavenumber = math.sqrt(epsb) * w_0 / c
		alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
		alpha = (alphasp**-1 - 1j*(2./3.)*wavenumber**3)**-1
		msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
		tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
		gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
		#print gamma_ret
		gamma_eV = gamma_ret*hbar/elec
		wsp = math.sqrt((wsp_0)**2); # sp frequency (rad/s) corrected for radiation damping
		print wsp
		for n in range (0,2*numPart):
			for m in range (0,2*numPart):
				if m == n: #if m and n are the same, the hammy gets the plasmon energy
					H[n,m] = 1
					#print H[n,m]
				elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				else: # all other dipoles are fair game
					R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
					p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
					p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = alpha*Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = alpha*(wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = alpha*(wavenumber**2)/(Rmag)
					space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
					space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
					ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
					gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
					H[n,m] = -np.real(ge) #this has the minus sign we need.
					#print H[n,m]
		w,v = scipy.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
		#print w
		idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
		eigenValues = w[idx] # sorting
		eigenVectors = v[:,idx] # sorting
		eigen=np.sqrt(eigenValues)*wsp*hbar/elec # the eigenvalues have units of energy^2, so we take the square root
		print eigen
	print eigen[2*numPart-1]
	numRings = 2
	E_ring = eigen[2*numPart-1]
	Loc_mag = [np.array([-e2,0,0]),np.array([e2,0,0])]
	near_vec = [np.array([2*e2,0]),np.array([e2,3*e1]),np.array([-e2,3*e1]),np.array([-2*e2,0]),np.array([-e2,-3*e1]),np.array([e2,-3*e1])]
	#next_near = [near_vec[0]+near_vec[1],near_vec[1]+near_vec[2],
	             #near_vec[2]+near_vec[3],near_vec[3]+near_vec[4],
	             #near_vec[4]+near_vec[5],near_vec[5]+near_vec[0]]
   	eigen_mag = np.zeros(numRings)
	alpha_m = 6*alphasp*0.5
	H_mag = np.zeros((numRings,numRings))
	E_plus = np.zeros((80,80))
	E_minus = np.zeros((80,80))
	tnn = nearest[r-1]/(2*E_ring)
	#tnnn = -(nn_sec[index*10]-ns_sec[index*10])*E_ring
	k_dot_near = 0
	k_dot_sec = 0
	kx = 1/(rij)
	ky = 1/(rij)
	k_scale = np.linspace(-math.pi,math.pi,80)
	for iii in range(0,80):
		for jjj in range(0,80):
			k_dot_near = 0
			k_dot_sec = 0
			k_tot = [k_scale[iii]*kx,k_scale[jjj]*ky]
			#for vec in range(0,3):
				#k_dot_sec = k_dot_sec + np.cos(np.dot(k_tot,next_near[vec]))
			for vec in near_vec:
				k_dot_near += (np.exp(1j*np.dot(k_tot,vec)))
				#print vec
				#print k_dot_near
			E_plus[iii,jjj] = np.sqrt(1+(tnn*abs(k_dot_near)))
			E_minus[iii,jjj] = np.sqrt(1-(tnn*abs(k_dot_near)))
	#np.savetxt('_'.join([str(r),'nm_E_plus.txt']),E_plus)
	#np.savetxt('_'.join([str(r),'nm_E_minus.txt']),E_minus)
	[iii, jjj] = np.meshgrid(np.linspace(-math.pi,math.pi,80),np.linspace(-math.pi,math.pi,80))
	fig = plt.figure()
	ax = fig.gca(projection = '3d')
	ax.plot_surface(iii,jjj,E_plus,color='red')
	ax.plot_surface(iii,jjj,E_minus,color='blue')
	plt.show()




