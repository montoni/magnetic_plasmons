import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt

# make six particles

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
frequency = np.linspace(1,6,1001) # we will need this later

'''nn_near = np.loadtxt('../../dielectric_study/NN_eps_1.txt')
ns_near = np.loadtxt('../../dielectric_study/NS_eps_1.txt')
nn_sec = np.loadtxt('../../unfused_twomer/NN_corners')
ns_sec = np.loadtxt('../../unfused_twomer/NS_corners')'''

''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''

for r in range(1,30):
	elec = 1.60217662e-19 # regular coulombs
	numPart = 6 #number of particles
	a0 = r*10**-7 #sphere radius in cm
	mie_index = (r-1)*10
	index = r-1
	ll = 1
	''' now determine geometry.'''
	''' begin with electric dipoles on the six-member ring '''
	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
	m_sep = 2*e2
	Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),
		   np.array([0,-2*e1]),np.array([e2,-e1]),np.array([e2,e1])]
	pol_vectors = [np.array([-1,0]),np.array([-0.5,-math.sqrt(3)/2]),np.array([0.5,-math.sqrt(3)/2]),
		   		   np.array([1,0]),np.array([0.5,math.sqrt(3)/2]),np.array([-0.5,math.sqrt(3)/2])]
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
	gamma = 0.05*elec/(hbar*16)
	wplasma = Eplasma/hbar; # plasma frequency (rad/s)
	''' First make epsilon(omega)'''

	freq_eps = epsinf - ((9.1)**2)/(np.square(frequency)+(0.05j/16)*frequency)
	#print freq_eps
	freq_alpha = (a0**3 )* ll*(freq_eps - epsb)/((ll*(freq_eps + epsb)) + epsb)
	alpha_mlwa = freq_alpha/(1-((2.0/3.0)*1j*(frequency/c*elec/hbar)**3)*freq_alpha - ((frequency/c * elec/hbar)**2)*freq_alpha/a0)
	wsp_0 = (mie_omegas[mie_index])*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = wsp_0
	eigen = np.zeros(2*numPart)
	count = 1
	S = 0
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
	alpha_eff = 1/((1/alpha_mlwa)-S)
	plt.figure(1)
	plt.plot(frequency,np.real(alpha_eff),frequency,np.imag(alpha_eff))#,frequency,freq_alpha)
	plt.show()
