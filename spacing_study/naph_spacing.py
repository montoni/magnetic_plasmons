import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''

elec = 1.60217662e-19 # regular coulombs
numPart = 10; #number of particles
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
epsb = 1
NN = []
NS = []
for dist in range(10,51):
	print dist
	r = 30
	a0 = r*10**-7; #sphere radius in cm

	''' now determine geometry.'''

	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * 2*(a0+(dist*10**-7)) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * 2*(a0+(dist*10**-7)); #long side of 30-60-90 triangle
	Loc = [np.array([0, e1]),np.array([-e2, 2*e1]),np.array([-2*e2, e1]),np.array([-2*e2, -e1]),np.array([-e2, -2*e1]),np.array([0, -e1]),np.array([e2, -2*e1]),np.array([2*e2, -e1]),np.array([2*e2 , e1]),np.array([e2 , 2*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.

	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = np.zeros((2*numPart,2*numPart)) #initialize Hammy with zeros, twenty by twenty in this case.

	'''More constants'''

	count = 1
	wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	'''initialize w_0 and eigen'''
	w_0 = 3*elec/hbar
	eigen = np.ones(2*numPart)
	for mode in range(0,2):
		while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-(mode+1)])) > 0.00000001:
			if count == 1:
				wsp = wsp_0
				count = count + 1
				w_0 = 0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
			print count
			alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
			msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(wsp**2) # I pulled this from the beats paper
			gamma_eV = gamma_ret*hbar/elec
			wsp = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,2*numPart):
				for m in range (n,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = (hbar*wsp/elec)**2
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = 1/Rmag**3 #this is the 1/r^3 term (static)
						r_squared = (w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						ge = (ch**2)*(r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
						print ge
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
		    #w_old = w_0
		    #w_0 = eigen[2*numPart-1]
		            
		if abs(np.sum(eigenVectors[:,(2*numPart)-(mode+1)])) <= 10**-10:
			NN.append(eigen[(2*numPart)-(mode+1)])
		else:
			NS.append(eigen[(2*numPart)-(mode+1)])
	
dist = np.arange(10,51)
plt.plot(dist,NN,dist,NS,linewidth=3)	
plt.show()
np.savetxt('_'.join([str(r),'nm_NN.txt']),NN)
np.savetxt('_'.join([str(r),'nm_NS.txt']),NS)
