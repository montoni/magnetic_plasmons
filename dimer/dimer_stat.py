import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''


r = 15
elec = 1.60217662e-19 # regular coulombs
numPart = 2; #number of particles
a0 = 15e-7; #sphere radius in cm

''' now determine geometry.'''

# make unit vectors in centimeters.
for number in range(0,41):
	dist = float(number+1)/4
	sep = (2+dist)*a0
	Loc = [np.array([sep/2, 0]),np.array([-sep/2, 0])]

	Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

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
	wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	'''initialize w_0 and eigen'''
	eigen = np.ones(2*numPart)


	alphasp = (a0**3)*(3/(epsinf+2)); # polarizability (cm^3)
	msp = (ch**2)/(alphasp*((wsp_0)**2)); # sp mass (grams)
	tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
	gamma_eV = gamma*hbar/elec
	for n in range (0,2*numPart):
		for m in range (n,2*numPart):
			if m == n: #if m and n are the same, the hammy gets the plasmon energy
				H[n,m] = (hbar*wsp_0/elec)**2
			elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
				H[n,m] = 0
			else: # all other dipoles are fair game
				R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
				Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
				nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
				p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
				p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
				r_cubed = 1/Rmag**3 #this is the 1/r^3 term (static)
				ge = (ch**2)*(r_cubed) * (3*p_nn_p - p_dot_p) #this is p dot E
				gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
				H[n,m] = -(ge/msp)*((hbar/elec)**2) #this has the minus sign we need.
	diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
	Ht = np.matrix.transpose(H) # this is the transpose of H
	Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
	Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
	w,v = scipy.linalg.eigh(Hfull) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
	idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
	eigenValues = w[idx] # sorting
	eigenVectors = v[:,idx] # sorting
	eigen=np.sqrt(eigenValues) # the eigenvalues have units of energy^2, so we take the square root
	print eigen
    #w_old = w_0
    #w_0 = eigen[2*numPart-1]
            
	print eigen #shows us the eigenvalues from largest to smallest
	print eigenVectors[:,(2*numPart)-1] # shows the lowest-energy mode
	print eigenVectors[:,(2*numPart)-2] # shows the second mode
	print eigenVectors[:,(2*numPart)-3] # shows the third mode
	print eigenVectors[:,(2*numPart)-4] # shows the fourth mode
	omega_plot = np.arange(3.5,4.0,0.0005)
	Lorentz_final = np.zeros(len(omega_plot))
	Lorentz = np.zeros((len(omega_plot),2*numPart))
	for omega_peak in range(0,2*numPart):
	    for counter in range(0,len(omega_plot)):
	        Lorentz[counter,omega_peak] = (1./(2.*math.pi)) * (gamma_eV)/((omega_plot[counter]-eigen[omega_peak])**2+(gamma_eV/2)**2)

	#plt.plot(omega_plot,np.sqrt(np.square(Lorentz[:,19])))
	#plt.show()
	np.savetxt('_'.join([str(dist),'nm_eigen.txt']),eigen)
	np.savetxt('_'.join([str(dist),'nm_spec.txt']),Lorentz)
	np.savetxt('_'.join([str(dist),'nm_lowest.txt']),eigenVectors[:,(2*numPart)-1])
	np.savetxt('_'.join([str(dist),'nm_second.txt']),eigenVectors[:,(2*numPart)-2])
