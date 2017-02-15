import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt


ll = 1
real_n = []
imag_n = []
mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
dipole = [] # doubly degenerate
no_dipole = [] # doubly degenerate
all_N = [] #all north
out_N_in_S = [] # that weird "excited state"
alt_NS = [] # maximally out of phase magnets
center = np.zeros(7,dtype=object)
frequency = np.linspace(1,6,1001)
for r in range(1,30):
	elec = 1.60217662e-19 # regular coulombs
	numPart = 24; #number of particles
	a0 = r*10**-7; #sphere radius in cm
	index = (r*10)
	''' now determine geometry.'''
	rad = 2.2*a0
	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
	Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),np.array([0,-2*e1]),
		   np.array([e2,-e1]),np.array([e2,e1]),np.array([2*e2,2*e1]),np.array([3*e2,e1]),
		   np.array([3*e2,-e1]),np.array([2*e2,-2*e1]),np.array([2*e2,-4*e1]),np.array([e2,-5*e1]),
		   np.array([0,-4*e1]),np.array([-e2,-5*e1]),np.array([-2*e2,-4*e1]),np.array([-2*e2,-2*e1]),
		   np.array([-3*e2,-e1]),np.array([-3*e2,e1]),np.array([-2*e2,2*e1]),np.array([-2*e2,4*e1]),
		   np.array([-e2,5*e1]),np.array([0, 4*e1]),np.array([e2,5*e1]),np.array([2*e2,4*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.
	center[0] = (Loc[0]+Loc[3])*0.5
	center[1] = (Loc[6]+Loc[9])*0.5
	center[2] = (Loc[9]+Loc[12])*0.5
	center[3] = (Loc[12]+Loc[15])*0.5
	center[4] = (Loc[15]+Loc[18])*0.5
	center[5] = (Loc[18]+Loc[21])*0.5
	center[6] = (Loc[21]+Loc[6])*0.5
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
	#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);

	freq_eps = epsinf - ((9.1)**2)/(np.square(frequency)+(0.05j)*frequency)
	#print freq_eps
	freq_alpha = 4*math.pi*(a0**3 )* ll*(freq_eps - epsb)/((ll*(freq_eps + epsb)) + epsb)
	alpha_mlwa = freq_alpha/(1-((1.0/(6.0*math.pi))*1j*((frequency/c)*(elec/hbar))**3)*freq_alpha - ((frequency/c * elec/hbar)**2)*freq_alpha/(a0*4*math.pi))

	wsp_0 = mie_omegas[index]*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = 3*elec/hbar
	eigen = np.ones(2*numPart)
	count = 1
	#for mode in range(0,7):
	while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-1])) > 0.00000001:
		if count == 1:
			w_0 = 0
			wsp = wsp_0
			count = count + 1
		else:
			w_0 = eigen[2*numPart-1]*elec/hbar
			count = count + 1
		print count
		print eigen[2*numPart-1]
		alphasp = (a0**3)*(3/(epsinf+2)); # polarizability (cm^3)
		msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
		tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
		gamma_ret = gamma+tau*(wsp**2) # I pulled this from the beats paper
		gamma_eV = gamma_ret*hbar/elec
		wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
		for n in range (0,2*numPart):
			for m in range (n,2*numPart):
				if m == n: #if m and n are the same, the hammy gets the plasmon energy
					H[n,m] = (hbar*wsp/elec)
				elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				else: # all other dipoles are fair game
					R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
					p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
					p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = alphasp/Rmag**3 #this is the 1/r^3 term (static)
					r_squared = (alphasp*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (alphasp*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					#space_exp = np.exp(1j*w_0*Rmag/c)
					space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
					space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
					ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
					gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
					H[n,m] = -(ge)*((hbar/elec)*wsp) #this has the minus sign we need.
		diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
		Ht = np.matrix.transpose(H) # this is the transpose of H
		Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
		Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
		w,v = scipy.linalg.eigh(Hfull) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
		idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
		eigenValues = w[idx] # sorting
		eigenVectors = v[:,idx] # sorting
		eigen=(eigenValues) # the eigenvalues have units of energy^2, so we take the square root
		#print eigen
	vec1 = np.reshape(eigenVectors[:,2*numPart-1],(numPart,2))
	vec2 = np.reshape(eigenVectors[:,2*numPart-2],(numPart,2))
	vec3 = np.reshape(eigenVectors[:,2*numPart-3],(numPart,2))
	vec4 = np.reshape(eigenVectors[:,2*numPart-4],(numPart,2))
	vec5 = np.reshape(eigenVectors[:,2*numPart-5],(numPart,2))
	vec6 = np.reshape(eigenVectors[:,2*numPart-6],(numPart,2))
	vec7 = np.reshape(eigenVectors[:,2*numPart-7],(numPart,2))

	for num_mode in [1,2,3,4,5,6,7]:
		vec = np.reshape(eigenVectors[:,2*numPart-num_mode],(numPart,2))
		mag_dipole = []
		if abs(np.sum(eigenVectors[:,2*numPart-num_mode])) < 1e-10:
			for cent in range(0,7):
				cr_pr = []
				for num in range(0,numPart):
					if num == numPart-1:
						idx = num - 1
					else:
						idx = num + 1
					if abs(np.linalg.norm(center[cent]-Loc[num]) - np.linalg.norm(Loc[num]-Loc[idx])) < 10**-10:
						cr_pr.append(np.cross(Loc[num]-center[cent],vec[num]))
					else:
						continue
				mag_dipole.append(np.sum(cr_pr))
			if abs(np.sum(mag_dipole)) < 1e-10:
				if abs(mag_dipole[0]) < 1e-10:
					if np.dot(vec[12],vec[21]) > 0:
						vec_alt_mag = vec
					else:
						pass
				else:
					pass
			else:
				if mag_dipole[0] * mag_dipole[1] > 0:
					vec_all_mag = vec
				else:
					vec_ring_mag = vec
		else:
			vec_ele = vec

		km = (elec/hbar)*frequency/c
	
	S = 0
	for n in range (0,numPart):
		for m in range (n,numPart):
			if n == m:
				S = S
		#elif m == n+1 and n%2 == 0:
		#	S = S
			else:
				R = Loc[(n)]-Loc[(m)] #pick out the location of each dipole, compute the distance between them
				Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
				nhat = (Loc[(n)]-Loc[(m)])/float(Rmag) #compute unit vector between dipoles
				p_dot_p = np.dot(vec_ring_mag[n],vec_ring_mag[m]) # this is one dipole dotted into the other
				p_nn_p = np.dot(vec_ring_mag[n],nhat)*np.dot(nhat,vec_ring_mag[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
				r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
				r_squared = 1j*(frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = np.square(frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
				exponent = np.exp((1j*Rmag*frequency*elec/hbar)/(c))
				S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)
	
	'''for n in range(0,numPart):
		for m in range(0,numPart):
			if m == n:
				S = S
			else:
				S = S + (3 + np.cos(2*math.pi*(m-n)/numPart))/((abs(np.sin(math.pi*(m-n)/numPart)))**3)'''

	km = (frequency*elec/hbar)/c

	alpha_eff_mag = ((1/alpha_mlwa) - S)**-1
	alpha_eff_mag = ((numPart*(km*rad)**2)/4)*(alpha_mlwa/(1-(alpha_mlwa*S/(4*math.pi)))) #*1j
	#alpha_eff_mag = ((4*epsb)/(numPart*(km*rad)**2)*(alpha_mlwa**-1 - S))**-1

	#alpha_eff = (4*epsb)/(numPart*(km**2)*(rad**2)*alpha_mlwa) - 1j*((km**3)/(6*math.pi) - (2*km)/(3*math.pi*numPart*rad**2)) + S/(16*math.pi*numPart*(km**2)*(rad**5))

	S = 0
	for n in range (0,numPart):
		for m in range (1,numPart):
			if n == m:
				S = S
			#elif m == n+1 and n%2 == 0:
			#	S = S
			else:
				R = Loc[(n)]-Loc[(m)] #pick out the location of each dipole, compute the distance between them
				Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
				nhat = (Loc[(n)]-Loc[(m)])/float(Rmag) #compute unit vector between dipoles
				p_dot_p = np.dot(vec_ele[n],vec_ele[m]) # this is one dipole dotted into the other
				p_nn_p = np.dot(vec_ele[n],nhat)*np.dot(nhat,vec_ele[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
				r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
				r_squared = 1j*(frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = np.square(frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
				exponent = np.exp((1j*Rmag*frequency*elec/hbar)/(c))
				S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)

	alpha_eff_ele = ((1/alpha_mlwa) - S)**-1
	alpha_eff_ele = ((4*epsb)/(numPart*(km*rad)**2)*(alpha_mlwa**-1 - S))**-1
	alpha_eff_ele = (1)*(alpha_mlwa/(1-(alpha_mlwa*S/(4*math.pi))))
	c_abs = 4*math.pi*(frequency/c)*np.imag(alpha_eff_mag)
	#print (2*math.pi*hbar*c)/(elec*np.power(np.imag(S[400]),(1./3.)))
	#plt.figure(1)
	#plt.plot(frequency,np.real(alpha_eff_mag),frequency,np.imag(alpha_eff_mag))#,frequency,freq_alpha)
	#plt.plot(frequency,np.real(alpha_eff_ele),frequency,np.imag(alpha_eff_ele))
	#plt.show()
	'''plt.figure(2)
	plt.plot(frequency,np.real(S),frequency,np.imag(S))
	plt.show()'''
	'''plt.figure()
	plt.plot(frequency,c_abs)
	plt.show()'''

	half_max = max(c_abs) / 2.
	#find when function crosses line half_max (when sign of diff flips)
	#take the 'derivative' of signum(half_max - Y[])
	d = np.sign(half_max - np.array(c_abs[0:-1])) - np.sign(half_max - np.array(c_abs[1:]))
	#plot(X,d) #if you are interested
	#find the left and right most indexes
	left_idx = np.where(d > 0)[0]
	right_idx = np.where(d < 0)[-1]
	fwhm = frequency[right_idx] - frequency[left_idx] #return the difference (full width)
	print fwhm
	volume = 6*a0**3
	#eps_cluster = (-volume - 2*alpha_eff)/(alpha_eff - volume)
	#print np.real(alpha_eff)
	res = zip(np.where(np.imag(alpha_eff_mag) == np.imag(alpha_eff_mag).max()))
	#print res
	#print frequency[res]

	### Okay! Time to calculate permeability! ###

	
	ring_density = 1/(2*a0*math.pi*rad**2)
	

	perm_eff = 1 + ((ring_density**-1)*(alpha_eff_mag**-1 + (1./(math.pi*6.))*1j*(km**3)))**-1
	eps_eff = 1 + ((ring_density**-1)*(alpha_eff_ele**-1 + (1./(math.pi*6.))*1j*(km**3)))**-1
	# - 1j/(km**3)
	#plt.figure()
	#plt.plot(frequency, np.real(perm_eff), frequency, np.imag(perm_eff))
	#plt.show()


	first_n = np.sqrt(eps_eff)
	second_n = np.sqrt(perm_eff)
	total_n = first_n * second_n
	plt.figure()
	plt.plot(frequency,np.real(total_n))
	plt.show()
	real_n.append(np.real(total_n))
	imag_n.append(np.imag(total_n))


	np.savetxt('real_n_sixmer_ring_mag',real_n)
	np.savetxt('imag_n_sixmer_ring_mag', imag_n)

	





