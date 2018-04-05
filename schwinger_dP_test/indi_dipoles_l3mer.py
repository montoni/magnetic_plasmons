import math
import numpy as np
import matplotlib.pyplot as plt

# these are unit vectors in each direction
zhat = [0,0,1]
yhat = [0,1,0]
xhat = [1,0,0]

p = 1
# flux screen distance
screen = 100000
#energies = np.array([3.61,3.6,3.56,3.52,3.42,3.3,3.17])
energies = np.array([3.25,1,1,3,2.8,1,1])
# some pre-defined stuff for loading data
#r_and_eV = [['1',3.61],['5',3.6],['10',3.56],['15',3.52],['20',3.42],['25',3.3],['30',3.17]]
r_and_eV = [['1',3.25],[],[],['15',3],['20',2.8]]

#NS_eig = np.loadtxt('../varying_particle_number/NS_eig.txt')
#NN_eig = np.loadtxt('../varying_particle_number/NN_eig.txt')


mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
#frequency = np.linspace(2,4,201) # we will need this later
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs

numPart = 14 #number of particles

me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; 
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
epsb = 1.77
frequency = np.linspace(2,4,201)
km = (frequency*elec/hbar)/c
NSN_vec = np.loadtxt('../chains/mode_vec_0.txt')
N_S_vec = np.loadtxt('../chains/mode_vec_1.txt')
NNN_vec = np.loadtxt('../chains/mode_vec_2.txt')
#dip_vec = np.loadtxt('../chains')
#vectors = [[vec_1],[vec_2],[vec_3]]
for sep in np.linspace(1,1,1):
	for r in [15]:
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print r
		#raw_input()
		# make unit vectors in centimeters.
		rij = (2+sep)*a0
		e1 = float(1)/float(2) * (rij) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (rij); #long side of 30-60-90 triangle
		Loc = [np.array([-3*e2,2*e1]),np.array([-4*e2,e1]),np.array([-4*e2,-e1]),np.array([-3*e2,-2*e1]),
			   np.array([0, e1]),np.array([-e2, 2*e1]),np.array([-2*e2, e1]),np.array([-2*e2, -e1]),
		       np.array([-e2, -2*e1]),np.array([0, -e1]),np.array([e2, -2*e1]),np.array([2*e2, -e1]),
		       np.array([2*e2 , e1]),np.array([e2 , 2*e1])]
		centers = [np.array([-3*e2,0]),np.array([-e2, 0]),np.array([e2,0])]
		places = Loc
		a0 = r*10**-7; #sphere radius in cm
		indind = int(r/5)
		omega_sp = mie_omegas[(r-1)*10] * elec/hbar
		wavenumber = energies[indind]*elec/(hbar*c)
		w_0 = energies[indind]*elec/hbar
		#3.42
		#print frequency
		freq_eps = epsinf - ((9.1)**2)/(np.square(frequency)+(0.05j)*frequency)
		#print freq_eps
		freq_alpha = (a0**3 )* (freq_eps - epsb)/(freq_eps + 2*epsb)
		alpha_mlwa = freq_alpha/(1-((2.0/3.0)*1j*((frequency/c)*(elec/hbar))**3)*freq_alpha - ((frequency/c * elec/hbar)**2)*freq_alpha/(a0))

		''' Load up some data '''

		data = np.loadtxt('../separation_polar/poyntyz_anth_'+str(r)+'_sep_2.6')
		#data = np.loadtxt('../water_dp_domega_3mers/water_yz_anth_'+str(r)+'_ene_'+str(r_and_eV[indind][1]))
		#data = np.loadtxt('water_spacing_yz_anth_15_ene_3')
		#data = np.loadtxt('poyntxz_rad_'+str(r)+'_ene_'+str(r_and_eV[indind][1]))
		#data = np.reshape(data,(len(data)/2,2),order='F')
		#data = np.sqrt(data[:,0]**2 + data[:,1]**2)
		data = data/data.max()
		
		alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
		alpha_ret = ((alphasp**-1) - 1j*(2./3.)*wavenumber**3)**-1
		#omega_sp = mie_omegas[(r-1)*10] * elec/hbar
		
		index = (r-1)*10
		''' now determine geometry.'''
		print index
		# make unit vectors in centimeters.
		rij = 3*a0
		rad = rij
		#print centers
		
		'''plt.figure()
		plt.scatter(xloc,yloc)
		plt.show()
		raw_input()'''
		#NN_vec = np.real(np.multiply(NN_vec,NN_p_dot_E))
		#NS_vec = np.real(np.multiply(NS_vec,NS_p_dot_E))
		# turn them into dipole moments
		#NS_vec = np.multiply(dipole_factor,NS_vec)
		#NN_vec = np.multiply(dipole_factor,NN_vec)
		#dip_vec = np.multiply(dipole_factor,dip_vec)

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
					p_dot_p = np.dot(NNN_vec[n],NNN_vec[m]) # this is one dipole dotted into the other
					p_nn_p = np.dot(NNN_vec[n],nhat)*np.dot(nhat,NNN_vec[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = 1j*(np.sqrt(epsb)*frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = np.square(np.sqrt(epsb)*frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					exponent = np.exp((1j*Rmag*np.sqrt(epsb)*frequency*elec/hbar)/(c))
					S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)
		#alpha_eff_mag = ((1/alpha_ret) - S)**-1
		alpha_eff_mag = ((numPart*(3*wavenumber*rad)**2)/4)*((alpha_mlwa**-1)-S)**-1
		print alpha_eff_mag
		raw_input()
		#maximum = alpha_eff_mag[np.where(alpha_eff_mag == alpha_eff_mag.max)]
		#print frequency[np.where(alpha_eff_mag == alpha_eff_mag.max())]
		#raw_input()
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
					p_dot_p = np.dot(N_S_vec[n],N_S_vec[m]) # this is one dipole dotted into the other
					p_nn_p = np.dot(N_S_vec[n],nhat)*np.dot(nhat,N_S_vec[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = 1j*(np.sqrt(epsb)*frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = np.square(np.sqrt(epsb)*frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					exponent = np.exp((1j*Rmag*np.sqrt(epsb)*frequency*elec/hbar)/(c))
					S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)
		#alpha_eff_ele = ((1/alpha_ret) - S)**-1
		alpha_eff_ele = (numPart)*(alpha_mlwa/(1-(alpha_mlwa*S)))
		

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
					p_dot_p = np.dot(NSN_vec[n],NSN_vec[m]) # this is one dipole dotted into the other
					p_nn_p = np.dot(NSN_vec[n],nhat)*np.dot(nhat,NSN_vec[m]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
					r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
					r_squared = 1j*(np.sqrt(epsb)*frequency*elec/hbar)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = np.square(np.sqrt(epsb)*frequency*elec/hbar)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
					exponent = np.exp((1j*Rmag*np.sqrt(epsb)*frequency*elec/hbar)/(c))
					S = S + ((r_unit * (p_dot_p - p_nn_p) + ((r_cubed - r_squared) * (3*p_nn_p - p_dot_p))) * exponent)
		alpha_eff_dip = (numPart)*(alpha_mlwa/(1-(alpha_mlwa*S)))
		print np.where(alpha_eff_mag == alpha_eff_mag.max())
		raw_input()
		dipole_factor_mag = hbar*alpha_eff_mag[np.where(alpha_eff_mag == alpha_eff_mag.max())]*omega_sp
		dipole_factor_ele = hbar*alpha_eff_ele[np.where(alpha_eff_mag == alpha_eff_mag.max())]*omega_sp
		dipole_factor_dip = hbar*alpha_eff_dip[np.where(alpha_eff_mag == alpha_eff_mag.max())]*omega_sp
		#dipole_factor_mag = hbar*alpha_eff_mag[143]*omega_sp
		#dipole_factor_ele = hbar*alpha_eff_ele[143]*omega_sp
		#dipole_factor_dip = hbar*alpha_eff_dip[143]*omega_sp
		#print alpha_eff_mag[np.where(alpha_eff_mag == alpha_eff_mag.max())]
		#print alpha_eff_ele[np.where(alpha_eff_mag == alpha_eff_mag.max())]
		#print alpha_eff_dip[np.where(alpha_eff_mag == alpha_eff_mag.max())]

		N_S_vec = (np.multiply(dipole_factor_ele,N_S_vec))
		NNN_vec = (np.multiply(dipole_factor_mag,NNN_vec))
		NSN_vec = (np.multiply(dipole_factor_dip,NSN_vec))
		
		p_m = 0#(np.multiply(p,dipole_factor_mag))
		p_e = (np.multiply(p,dipole_factor_ele))

		#for x in range(0,numPart):
			#electric = np.multiply([0,1,0],np.exp(1j*wavenumber*np.sqrt(places[x][0]**2 + places[x][1]**2)))
			#NS_p_dot_E += np.dot(NS_vec[x],electric)
			#NN_p_dot_E += np.dot(NN_vec[x],electric)
		d_list = [[],[]]
		mu_list = [[],[]]
		# "drive" dipole moments with a "field"
		
		for cent in [0,1]:
			origin = centers[cent]
			mag_dipole = 0
			elec_dipole = []
			for num in range(0,numPart):
				#print np.linalg.norm(origin - places[num])
				#raw_input()
				if np.linalg.norm(origin - places[num]) > 3.1*a0:
					#print num
					#print "passed"
					pass
				else:
					#print num
					#print "not passed"
					elec_dipole.append(N_S_vec[num])
					#elec_dipole.append(dip_vec[num])
					mag_dipole += np.cross( origin - places[num],NNN_vec[num])
			#raw_input()
			d_list[cent].append(np.sum(elec_dipole))
			mu_list[cent].append(mag_dipole*energies[indind]*elec/(hbar*2*c))
		print d_list
		print mu_list
		'''m_dot_B = 0+0*1j
		p_dot_E = 0+0*1j
		for x in [0,1]:
			magnetic = np.exp(1j*wavenumber*np.linalg.norm(centers[x]))
			#print mu_list[x]
			m_dot_B += np.multiply(mu_list[x],magnetic)
			p_dot_E += np.multiply(d_list[x],magnetic)
		d_list = (np.multiply(p_dot_E,d_list))
		mu_list = (np.multiply(m_dot_B,mu_list))
		print d_list
		print mu_list'''
		dP_dOmega = np.zeros((37,1))
		#centers = [np.array([-cent_dist,0,0]),np.array([cent_dist,0,0])]

		theta_count = 0
		for theta in np.linspace(0,2*math.pi,37):
			nhat = [0*np.cos(theta),np.cos(theta),np.sin(theta)]
			
			electric_dipoles = ((abs(p_e))**2)*(1-np.dot(nhat,yhat)**2)
			magnetic_dipoles = ((abs(p_m))**2)*(1-np.dot(nhat,zhat)**2)*(wavenumber**2)*(3*(rad)**2)
			interference_1 = (p_e*np.conj(p_m)+p_m*np.conj(p_e))*wavenumber*(3*rad)*np.dot(nhat,xhat)
			dP_dOmega[theta_count] += np.real((electric_dipoles + magnetic_dipoles + interference_1))
			theta_count += 1
		#print dP_dOmega
		#raw_input()
		theta = np.linspace(0,2*math.pi,37)
		theta2 = np.linspace(0,2*math.pi,361)
		plt.figure()
		plt.polar(theta,dP_dOmega/dP_dOmega.max(),marker='o',linewidth=0)
		plt.polar(theta2,data)
		plt.legend(['theory','sim'])
		plt.savefig('dp_do_vac_spaced_yz_l3mer_'+str(r)+'.pdf')
		plt.show()