import numpy as np
import scipy
import math
import matplotlib.pyplot as plt
import string

# make six particles
for den in np.linspace(0.4,0.4,1):
	
	frequency = np.linspace(0,2,201) # we will need this later
	real_n = []
	imag_n = []

	''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''

	for r in np.linspace(1.5,10,35):
		dielectrics = np.loadtxt('density_' + str(den) + '/epsilon_den_' + str(den) + '_rad_' + str(r))
		#print dielectrics[2]
		freq_index, maximum = max(enumerate(dielectrics[2]), key=lambda item: item[1])
		print dielectrics[0][freq_index]
		elec = 1.60217662e-19 # regular coulombs
		numPart = 6 #number of particles
		a0 = r*10**-7 #sphere radius in cm
		ll = 1
		''' now determine geometry.'''
		''' begin with electric dipoles on the six-member ring '''
		# make unit vectors in centimeters.
		e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
		mag_x = 2*e2
		mag_y = 3*e1
		rad = 2.2*a0
		Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),np.array([0,-2*e1]),np.array([e2,-e1]),np.array([e2,e1])]
		#pol_vectors = [np.array([-1/math.sqrt(3),0]),np.array([-0.5/math.sqrt(3),-math.sqrt(3)/2/math.sqrt(3)]),np.array([0.5/math.sqrt(3),-math.sqrt(3)/2/math.sqrt(3)]),
			   		   #np.array([1/math.sqrt(3),0]),np.array([0.5/math.sqrt(3),math.sqrt(3)/2/math.sqrt(3)]),np.array([-0.5/math.sqrt(3),math.sqrt(3)/2/math.sqrt(3)])]
		vec_mag = [np.array([-1,0]),np.array([-0.5,-math.sqrt(3)/2]),np.array([0.5,-math.sqrt(3)/2]),
			   		   np.array([1,0]),np.array([0.5,math.sqrt(3)/2]),np.array([-0.5,math.sqrt(3)/2])]
		vec_ele = [[0,1],[0,1],[0,1],[0,1],[0,1],[0,1]]
		Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

		'''More constants'''
		
		me = 9.10938291e-28; # electron mass in g
		ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
		hbar = 1.054571726e-34; # modified Planck in J*s
		c = 2.99e10; # speed of light in cm/s
		eps0 = 8.85418782e-12; # permittivity of free space
		epsb = 1; # background dielectric - right now we're in vacuum
		epsinf = 3.71; # does this have units?
		'''Properties for silver.'''
		Eplasma = 0.928317; # eV
		gamma = (elec/hbar)*0.075+(0.1*0.553)/(a0*1e-2) #Hz
		print gamma
		wplasma = Eplasma/hbar; # plasma frequency (rad/s)
		''' First make epsilon(omega)'''
		km = (frequency/c)*(elec/hbar)
		freq_eps = dielectrics[1] + 1j*dielectrics[2]
		#print freq_eps
		freq_alpha = 4*math.pi*(a0**3 )* ll*(freq_eps - epsb)/((ll*(freq_eps + epsb)) + epsb)
		alpha_mlwa = freq_alpha/(1-((2.0/3.0)*1j*((frequency/c)*(elec/hbar))**3)*freq_alpha )
		alpha_rad = (freq_alpha**-1 - 1j*(1./(math.pi*6.))*km**3 - (km**2)/(4*math.pi*a0))**-1
		wsp_0 = (dielectrics[0][freq_index])*elec/hbar
		print wsp_0
		#initialize w_0 and eigen
		w_0 = wsp_0
		alphasp_stat = ((a0**3)*3)/(epsinf+2)
		eigen = np.zeros(2*numPart)
		count = 1
		H = np.zeros((2*numPart,2*numPart))
		while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-1])) > 0.00000001:
			if count == 1:
				w_0 = wsp_0
				count = count + 1
				wsp = wsp_0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-1]*elec/hbar
				#if w_0 > wsp:
				#	w_0 = 3.5*elec/hbar
				msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
				tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
				gamma_ret = gamma+tau*(w_0**2) 
				#print gamma_ret
				gamma_eV = gamma_ret*hbar/elec
				wsp = math.sqrt((wsp_0)**2 - (gamma/2)**2);
				print wsp
			print eigen[2*numPart-1]
			eps_sp = dielectrics[1][freq_index] + 1j*dielectrics[2][freq_index]
			#print eps_sp
			alphasp = (a0**3)*(eps_sp - epsb)/(eps_sp + 2*epsb); # polarizability (cm^3)
			alphasp_mlwa = alphasp/(1 - (2./3.)*1j*(w_0/c)*alphasp-(alphasp/a0)*(w_0/c)**2)
			print wsp*hbar/elec
			for n in range (0,2*numPart):
				for m in range (n,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = (hbar*wsp/elec)
						#print H[n,m]
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alphasp_stat*(Rmag**-3) #this is the 1/r^3 term (static)
						r_squared = (alphasp_stat*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alphasp_stat*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -(ge*wsp)*((hbar/elec)) #this has the minus sign we need.
						#print H[n,m]
			diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
			Ht = np.matrix.transpose(H) # this is the transpose of H
			Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			print Hfull
			w,v = np.linalg.eigh(Hfull) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			#print w
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting
			
			eigen=(eigenValues) # redefine
			lowest = eigen[2*numPart-1]
		print eigen
		print eigenVectors
		vec_mag = np.reshape(eigenVectors[:,(2*numPart)-1],(numPart,2))
		vec_ele = np.reshape(eigenVectors[:,(2*numPart)-2],(numPart,2))
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
					p_dot_p = np.dot(vec_mag[n],vec_mag[m]) # this is one dipole dotted into the other
					p_nn_p = np.dot(vec_mag[n],nhat)*np.dot(nhat,vec_mag[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
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

		#alpha_eff_mag = ((1/alpha_mlwa) - S)**-1
		#print (km*rad)**2
		alpha_eff_mag = ((numPart*(km*rad)**2)/4)*(alpha_rad/(1-(alpha_rad*S/(4*math.pi)))) #*1j
		#alpha_eff_mag = ((4*epsb)/(numPart*(km*rad)**2)*(alpha_mlwa**-1 - S))**-1 *(km*rad)**2

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

		#alpha_eff_ele = ((1/alpha_mlwa) - S)**-1
		#alpha_eff_ele = ((4*epsb)/(numPart*(km*rad)**2)*(alpha_mlwa**-1 - S))**-1
		alpha_eff_ele = (numPart)*(alpha_rad/(1-(alpha_rad*S/(4*math.pi))))
		#c_abs = 4*math.pi*(frequency/c)*np.imag(alpha_eff_mag)
		#print (2*math.pi*hbar*c)/(elec*np.power(np.imag(S[400]),(1./3.)))
		plt.figure(1)
		plt.plot(frequency,np.real(alpha_eff_mag),frequency,np.imag(alpha_eff_mag))
		plt.plot(frequency,np.real(alpha_eff_ele),frequency,np.imag(alpha_eff_ele))
		plt.show()
		'''plt.figure(2)
		plt.plot(frequency,np.real(S),frequency,np.imag(S))
		plt.show()'''
		'''plt.figure()
		plt.plot(frequency,c_abs)
		plt.show()'''

		### Okay! Time to calculate permeability! ###

		
		ring_density = 1/(2*a0*math.pi*rad**2)
		

		perm_eff = 1 + alpha_eff_mag*((ring_density**-1)*(1 + (1./(math.pi*6.))*1j*(alpha_eff_mag*km**3)))**-1
		eps_eff = 1 + ((ring_density**-1)*(alpha_eff_ele**-1 + (1./(math.pi*6.))*1j*(km**3)))**-1
		# - 1j/(km**3)
		plt.figure()
		plt.plot(frequency, np.real(perm_eff), frequency, np.imag(perm_eff))
		plt.plot(frequency, np.real(eps_eff), frequency, np.imag(eps_eff))
		plt.show()


		first_n = np.sqrt(eps_eff)
		second_n = np.sqrt(perm_eff)
		total_n = first_n * second_n
		plt.figure()
		plt.plot(frequency,np.real(total_n))
		plt.show()
		real_n.append(np.real(total_n))
		imag_n.append(np.imag(total_n))
		np.savetxt('density_'+str(den)+'/real_n_onemer',np.real(total_n))
		np.savetxt('density_'+str(den)+'/imag_n_onemer',np.imag(total_n))


		





