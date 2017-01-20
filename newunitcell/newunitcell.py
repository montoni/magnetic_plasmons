import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

ll = 1
mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
north_north = np.loadtxt('../dielectric_study/NN_right.txt')
north_south = np.loadtxt('../dielectric_study/NS_right.txt')
frequency = np.linspace(1,6,1001)
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
NN = []
NS = []
part_positions = []
for r in range(1,30):
	elec = 1.60217662e-19 # regular coulombs
	 #number of particles
	a0 = r*10**-7 #sphere radius in cm
	mie_index = (r-1)*10
	index = r-1
	flag = 'yes'
	''' now determine geometry.'''
	''' begin with electric dipole on the six-member ring '''
	# make unit vectors in centimeters.
	ey = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	ex = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
	mag_rad = np.sqrt(np.square(ex)+np.square(ey))
	uy = 6*ey
	ux = 4*ex
	unitcell = 1
	Loc = [np.array([-2*ex,-2*ey]),np.array([-ex,-ey]),np.array([0,-2*ey]),np.array([-ex,ey]),
	       np.array([0,2*ey]),np.array([ex,-ey]),np.array([ex,ey]),np.array([2*ex,2*ey])]

	part_positions = []
	numPart = unitcell*8
	for num in range(0,unitcell): 
		if flag == 'yes':
			for mun in range(0,unitcell):
				Loc = [np.array([-2*ex+num*ux,-2*ey+mun*uy]),
				   	   np.array([-ex+num*ux,-ey+mun*uy]),
				       np.array([num*ux,-2*ey+mun*uy]),
				       np.array([-ex+num*ux,ey+mun*uy]),
				       np.array([num*ux,2*ey+mun*uy]),
				       np.array([ex+num*ux,-ey+mun*uy]),
				       np.array([ex+num*ux,ey+mun*uy]),
				       np.array([2*ex+num*ux,2*ey+mun*uy])]
				
				for pos in range(0,np.shape(Loc)[0]):
					part_positions.append(Loc[pos]) 
				print part_positions
	Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction


	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = np.zeros((2*numPart,2*numPart)) #initialize Hammy with zeros

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
		alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
		msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
		tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
		gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
		#print gamma_ret
		gamma_eV = gamma_ret*hbar/elec
		wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
		print wsp
		for n in range (0,2*numPart):
			for m in range (n,2*numPart):
				if m == n: #if m and n are the same, the hammy gets the plasmon energy
					H[n,m] = (hbar*wsp/elec)
					#print H[n,m]
				elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				else: # all other dipoles are fair game
					print n
					print m
					R = part_positions[(n/2)]-part_positions[(m/2)] #pick out the location of each dipole, comute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					nhat = (part_positions[(n/2)]-part_positions[(m/2)])/float(Rmag) #compute unit vector between dipoles
					p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
					p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = alphasp*Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = (alphasp*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (alphasp*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					#space_exp = np.exp(1j*w_0*Rmag/c)
					space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
					space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
					ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
					gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
					H[n,m] = -(ge)*(wsp*(hbar/elec)) #this has the minus sign we need.
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
		eigen=(eigenValues) # the eigenvalues have units of energy^2, so we take the square root
		print eigen


	x,y = zip(*part_positions)
	vec = np.reshape(eigenVectors[:,(2*numPart)-1], (numPart,2))
	#mag_dipole = (np.sum(np.linalg.norm(vec,axis=1)) + np.linalg.norm(vec[:,0],axis=1) np.linalg.norm(vec[:,7],axis=1))*eigen[2*numPart-1]*mag_rad/(2*c)
	#numRings = 2
	#mag_loc = [np.array([0,0]),np.array(ux,0)]

	### now compute polarizability of a single ring ###
	freq_eps = epsinf - ((9.1)**2)/(np.square(frequency)+(0.05j/16)*frequency)
	#print freq_eps
	freq_alpha = (a0**3 )* ll*(freq_eps - epsb)/((ll*(freq_eps + epsb)) + epsb)
	alpha_mlwa = freq_alpha/(1-((2.0/3.0)*1j*(frequency/c*elec/hbar)**3)*freq_alpha - ((frequency/c * elec/hbar)**2)*freq_alpha/a0)
	
	'''initialize w_0 and eigen'''
	
	
	S = 0
	#u,v = zip(*vec)
	pol_vectors = vec
	norm_factor = 1./np.sqrt(np.square(pol_vectors[:,0]) + np.square(pol_vectors[:,1]))

	print pol_vectors
	print norm_factor
	for n in range (0,numPart):
		for m in range (n,numPart):
			if n == m:
				S = S
			#elif m == n+1 and n%2 == 0:
				#S = S
			else:
				R = Loc[(n)]-Loc[(m)] #pick out the location of each dipole, compute the distance between them
				Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
				nhat = (Loc[(n)]-Loc[(m)])/float(Rmag) #compute unit vector between dipoles
				p_dot_p = np.dot(pol_vectors[n],pol_vectors[m]) # this is one dipole dotted into the other
				p_nn_p = np.dot(pol_vectors[n],nhat)*np.dot(nhat,pol_vectors[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
				r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
				r_squared = 1j*(frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = np.square(frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
				exponent = np.exp((1j*Rmag*frequency*elec/hbar)/(c))
				S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)
	alpha_eff = 1/((1/alpha_mlwa)-S*2)
	plt.figure(1)
	plt.plot(frequency,np.real(alpha_eff),frequency,np.imag(alpha_eff))# frequency,freq_alpha)
	plt.plot(frequency,np.real(alpha_mlwa),frequency,np.imag(alpha_mlwa))
	plt.show()


	#print Loc[:]
	
	#y = Loc[:,1]
	#u,v = zip(*vec)
	#v = vec[:,1]
	#plt.title(' '.join(['radius =',str(r)]))
	#plt.quiver(x,y,u,v)
	#plt.show()
	