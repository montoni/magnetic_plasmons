import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
nn_near = np.loadtxt('../dielectric_study/NN_right.txt')
ns_near = np.loadtxt('../dielectric_study/NS_right.txt')
nn_sec = np.loadtxt('../unfused_twomer/NN_sec_near')
ns_sec = np.loadtxt('../unfused_twomer/NS_sec_near')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
NNN = []
NSN = []
N_S = []

for r in range(10,301):
	elec = 1.60217662e-19 # regular coulombs
	numPart = 6 #number of particles
	a0 = .1*r*10**-7 #sphere radius in cm
	mie_index = r-10
	index = r-10

	''' now determine geometry.'''
	''' begin with electric dipole on the six-member ring '''
	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
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
		#print eigen[2*numPart-1]
		alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
		msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
		tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
		gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
		#print gamma_ret
		gamma_eV = gamma_ret*hbar/elec
		wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
		#print wsp
		for n in range (0,2*numPart):
			for m in range (n,2*numPart):
				if m == n: #if m and n are the same, the hammy gets the plasmon energy
					H[n,m] = (hbar*wsp/elec)**2
					#print H[n,m]
				elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				else: # all other dipoles are fair game
					R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
					p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
					p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = (w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					#space_exp = np.exp(1j*w_0*Rmag/c)
					space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
					space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
					ge = (ch**2)*(r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
					gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
					H[n,m] = -(ge/msp)*((hbar/elec)**2) #this has the minus sign we need.
					#print H[n,m]
		diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
		Ht = np.matrix.transpose(H) # this is the transpose of H
		Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
		Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
		#print Hfull
		w,v = scipy.linalg.eigh(Hfull) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
		#print w
		idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
		eigenValues = w[idx] # sorting
		eigenVectors = v[:,idx] # sorting
		eigen=np.sqrt(eigenValues) # the eigenvalues have units of energy^2, so we take the square root
		#print eigen
	print eigen[2*numPart-1]
	#raw_input("press Enter to continue...")
	numRings = 3
	E_ring = eigen[2*numPart-1]
	omega_ring = E_ring*elec/hbar
	Loc_mag = [np.array([-2*e2,0,0]),np.array([0,0,0]),np.array([2*e2,0,0])]
   	eigen_mag = np.zeros(numRings)
	alpha_m = 6*alphasp*0.5
	H_mag = np.zeros((numRings,numRings))
	w_mag = omega_ring
	coupling_near = nn_near[index]-ns_near[index]
	coupling_sec = nn_sec[index]-ns_sec[index]
	for mode in range(1,numRings+1):
		count = 1
		while np.sqrt(np.square(w_mag*hbar/elec - eigen_mag[numRings-mode])) > 0.00000001:
			if count == 1:
				w_mag = 0*elec/hbar
				count = count + 1
				w_mag_0 = omega_ring
			else:
				count = count + 1
				w_mag = eigen_mag[numRings-mode]*elec/hbar
			print count
			m_mag = (ch**2)/(alpha_m*((w_mag_0)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*m_mag*c**3) # radiation damping time
			gamma_ret_mag = gamma+tau*(w_mag_0**2)
			#print gamma_ret
			gamma_eV = gamma_ret*hbar/elec
			w_mag_0 = np.sqrt((omega_ring)**2 - (gamma_ret_mag/2)**2);  #sp frequency (rad/s) corrected for radiation damping
			Q_mag = [0,0,(6*w_mag_0*2*e1)/(2*c)]
			#
			for n in range (0,numRings):
				for m in range (n,numRings):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H_mag[n,m] = (hbar*w_mag_0/elec)**2
						#print H[n,m]
					else: # all other dipoles are fair game
						R = Loc_mag[n]-Loc_mag[m] #pick out the location of each dipole, compute the distance between them
						Rmag = math.sqrt(R[0]**2 + R[1]**2)#math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						if 1.9*e2 < Rmag < 2.1*e2:
							H_mag[n,m] = ((coupling_near))*(w_mag_0*hbar/elec)
						if 3.8*e2 < Rmag < 4.2*e2:
							H_mag[n,m] = ((coupling_sec))*(w_mag_0*hbar/elec)
						#print H_mag
						#raw_input("press enter")
			diag_mag = np.diag(np.diag(H_mag)) # this produces a matrix of only the diagonal terms of H
			Ht_mag = np.matrix.transpose(H_mag) # this is the transpose of H
			Hedit_mag = diag_mag - Ht_mag # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			Hfull_mag = H_mag - Hedit_mag # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			#print Hfull_mag
			w_m,v_m = scipy.linalg.eigh(Hfull_mag) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			#print w
			idx = w_m.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigenValues_mag = w_m[idx] # sorting
			eigenVectors_mag = v_m[:,idx] # sorting
			eigen_mag=np.sqrt(eigenValues_mag) # the eigenvalues have units of energy^2, so we take the square root
			#print eigen_mag
			print eigenVectors_mag
		mag_vec = eigenVectors_mag[:,mode-1]
		#print mag_vec
		if mag_vec[0]*mag_vec[2] < 0:
			N_S.append(eigen_mag[mode-1])
		elif mag_vec[0]*mag_vec[1] < 0:
			NSN.append(eigen_mag[mode-1])
		else:
			NNN.append(eigen_mag[mode-1])

print len(NNN)
print len(NSN)
print len(N_S)
'''np.savetxt('NN_wrong.txt',NN)
np.savetxt('NS_wrong.txt',NS)'''

r = np.linspace(1,30,291)
plt.plot(r,NNN,r,NSN,r,N_S,linewidth=3)	
plt.legend(['NNN','NSN','N_S'])
plt.show()