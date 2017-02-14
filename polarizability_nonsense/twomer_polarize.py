import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import string

# make six particles

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
frequency = np.linspace(1,6,1001) # we will need this later
real_n = []
imag_n = []
'''nn_near = np.loadtxt('../../dielectric_study/NN_eps_1.txt')
ns_near = np.loadtxt('../../dielectric_study/NS_eps_1.txt')
nn_sec = np.loadtxt('../../unfused_twomer/NN_corners')
ns_sec = np.loadtxt('../../unfused_twomer/NS_corners')'''

''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''

for r in range(1,30):
	elec = 1.60217662e-19 # regular coulombs
	numPart = 10 #number of particles
	a0 = r*10**-7 #sphere radius in cm
	mie_index = (r*10)
	index = r-1
	ll = 1
	rad = 2.2*a0
	''' now determine geometry.'''
	''' begin with electric dipoles on the six-member ring '''
	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
	m_sep = 2*e2
	Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),
		   np.array([0,-2*e1]),np.array([e2,-e1]),np.array([e2,e1]),
		   np.array([2*e2,2*e1]),np.array([3*e2,e1]),np.array([4*e2,-e1]),
		   np.array([3*e2,-2*e1])]
	vec_mag = [[-1,0],[-0.5,-math.sqrt(3)/2],[0.5,-math.sqrt(3)/2],[1,0],[1,0],[-1,0],[-1,0],[-0.5,math.sqrt(3)/2],[0.5,math.sqrt(3)/2],[1,0]]
	vec_ele = [[-1,0],[-0.5,-math.sqrt(3)/2],[0.5,-math.sqrt(3)/2],[1,0],[0,1],[0,1],[1,0],[0.5,-math.sqrt(3)/2],[-0.5,-math.sqrt(3)/2],[-1,0]]
	Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

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
	gamma = 0.05*elec/(hbar)
	wplasma = Eplasma/hbar; # plasma frequency (rad/s)
	''' First make epsilon(omega)'''

	freq_eps = epsinf - ((9.1)**2)/(np.square(frequency)+(0.05j)*frequency)
	#print freq_eps
	freq_alpha = (a0**3 )* ll*(freq_eps - epsb)/((ll*(freq_eps + epsb)) + epsb)
	alpha_mlwa = freq_alpha/(1-((2.0/3.0)*1j*((frequency/c)*(elec/hbar))**3)*freq_alpha - ((frequency/c * elec/hbar)**2)*freq_alpha/a0)
	wsp_0 = (mie_omegas[mie_index])*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = wsp_0
	alphasp_stat = ((a0**3)*3)/(epsinf+2)

	km = (elec/hbar)*frequency/c
	
	S = 0
	for n in range (0,1):
		for m in range (n,numPart):
			if n == m:
				S = S
		#elif m == n+1 and n%2 == 0:
		#	S = S
			else:
				R = Loc[(n)]-Loc[(m)] #pick out the location of each dipole, compute the distance between them
				Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
				nhat = (Loc[(n)]-Loc[(m)])/float(Rmag) #compute unit vector between dipoles
				p_dot_p = np.dot(vec_mag[n],vec_mag[m]) # this is one dipole dotted into the other
				p_nn_p = np.dot(vec_mag[n],nhat)*np.dot(nhat,vec_mag[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
				r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
				r_squared = 1j*(frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = np.square(frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
				exponent = np.exp((1j*Rmag*frequency*elec/hbar)/(c))
				S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)

	alpha_eff_mag = ((1/alpha_mlwa) - S)**-1
	alpha_eff_mag = (numPart*km*rad/2)*(alpha_mlwa/(1-alpha_mlwa*S)) #*1j
	#alpha_eff_mag = ((4*epsb)/(numPart*(km*rad)**2)*(alpha_mlwa**-1 - S))**-1

	#alpha_eff = (4*epsb)/(numPart*(km**2)*(rad**2)*alpha_mlwa) - 1j*((km**3)/(6*math.pi) - (2*km)/(3*math.pi*numPart*rad**2)) + S/(16*math.pi*numPart*(km**2)*(rad**5))

	S = 0
	for n in range (0,1):
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
	alpha_eff_ele = (numPart*km*rad/2)*(alpha_mlwa/(1-alpha_mlwa*S))
	c_abs = 4*math.pi*(frequency/c)*np.imag(alpha_eff_mag)
	#print (2*math.pi*hbar*c)/(elec*np.power(np.imag(S[400]),(1./3.)))
	plt.figure(1)
	plt.plot(frequency,np.real(alpha_eff_mag),frequency,np.imag(alpha_eff_mag))#,frequency,freq_alpha)
	plt.plot(frequency,np.real(alpha_eff_ele),frequency,np.imag(alpha_eff_ele))
	plt.show()
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

	
	ring_density = 1/(2*2*a0*rad*e2/2)
	

	perm_eff = 1 + ((ring_density**-1)*((alpha_eff_mag**-1)) - (1./3.))**-1
	eps_eff = 1 + ((ring_density**-1)*((alpha_eff_ele**-1)) - (1./3.))**-1
	#plt.figure()
	#plt.plot(frequency, np.real(perm_eff), frequency, np.imag(perm_eff))
	#plt.show()


	first_n = np.sqrt(eps_eff)
	second_n = np.sqrt(perm_eff)
	total_n = first_n * second_n
	plt.figure()
	plt.plot(frequency,np.real(total_n),frequency,np.imag(total_n))
	plt.show()

	real_n.append(np.real(total_n))
	imag_n.append(np.imag(total_n))

	np.savetxt('real_n_twomer',real_n)
	np.savetxt('imag_n_twomer', imag_n)

	





