import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs

numPart = 10; #number of particles
omegas = []
me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; 
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar*16)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
epsb = 1
NN = []
NS = []
modes = [[],[]]
interaction = [[],[]]
r = 20
for epsb in np.linspace(1,3,21):
	a0 = r*10**-7; #sphere radius in cm
	alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
	index = r
	''' now determine geometry.'''
	print epsb
	# make unit vectors in centimeters.
	rij = 3*a0
	part_per_ring = numPart/2 + 1
	theta = 2*math.pi/part_per_ring
	phi = theta/2.
	#print theta
	#print phi
	#print np.cos(phi)
	#print np.sin(phi)
	cent_dist = rij/(2*np.tan(phi))
	part_to_cent = math.sqrt((rij/2)**2 + (cent_dist)**2)
	centers = [np.array([-cent_dist,0]),np.array([cent_dist,0])]
	#print centers
	places = []
	for num in range(part_per_ring-1):
		if part_per_ring%2 == 0:
			#print 'hello'
			places.append(centers[0] + np.array([(part_to_cent*np.cos(phi+theta*(num))),(part_to_cent*np.sin(phi+theta*(num)))]))
			places.append(centers[1] + np.array([(part_to_cent*np.cos(phi+theta*(num+np.ceil(part_per_ring/2.)))),(part_to_cent*np.sin(phi+theta*(num+np.ceil(part_per_ring/2.))))]))
		elif part_per_ring%2 == 1:
			#print 'good bye'
			places.append(centers[0] + np.array([(part_to_cent*np.cos(phi+theta*(num))),(part_to_cent*np.sin(phi+theta*(num)))]))
			places.append(centers[1] + np.array([(part_to_cent*np.cos(theta*(num+np.ceil(part_per_ring/2.)))),(part_to_cent*np.sin(theta*(num+np.ceil(part_per_ring/2.))))]))
	'''Loc = np.zeros((numPart,2))
	count = 0 
	for x in places:
		print x
		if np.isclose(x.all(),Loc[:].any()):
			pass
		else:
			Loc[count] = x'''


	#Loc = np.unique(Loc)
	#print places
	Loc = places
	xloc, yloc = zip(*Loc)
	#plt.scatter(xloc,yloc)
	#raw_input()
	#plt.show()
	#raw_input()
	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = (np.zeros((2*numPart,2*numPart),dtype=float))#initialize Hammy with zeros, twenty by twenty in this case.

	'''More constants'''

	count = 1
	#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	wsp_0 = (mie_omegas[index*10 - 10])*(math.sqrt(epsinf+2)/math.sqrt(epsinf+2*epsb))    *elec/hbar
	omegas.append(wsp_0*hbar/elec)
	'''initialize w_0 and eigen'''
	w_0 = 0
	eigen = np.ones(2*numPart)
	for mode in range(0,2):
		while np.sqrt(np.square(w_0*hbar/elec - eigen[(2*numPart)-(mode+1)])) > 0.00001:
			if count == 1:
				wsp = wsp_0
				count = count + 1
				w_0 = 0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
			wavenumber = (w_0*math.sqrt(epsb))/(c)
			alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
			msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
			gamma_eV = gamma_ret*hbar/elec
			#wsp = wsp_0
			#wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2) # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,2*numPart):
				for m in range (0,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1#(hbar/elec)*wsp
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, compute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*wavenumber)/(Rmag**2) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*wavenumber**2)/(Rmag) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(wavenumber*Rmag) #this is the real part of the e^ikr
						space_sin = np.sin(wavenumber*Rmag) #this is the imaginary part of the e^ikr
						exponent = np.exp(1j*wavenumber*Rmag)
						ge = ((r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))) * exponent #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -(ge)/epsb#*(hbar/elec)*wsp #this has the minus sign we need.
						H[m,n] = np.conj(-ge)/epsb
						

					

			'''diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
			Ht = np.matrix.transpose(H) # this is the transpose of H
			Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			#print Hfull
			#raw_input()'''
			w,v = scipy.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			#print idx
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting

			eigen=((hbar/elec)*wsp)*(np.sqrt(eigenValues))# the eigenvalues have units of energy^2, so we take the square root
		vec = np.reshape(eigenVectors[:,2*numPart - (mode+1)],[numPart,2])
		#print vec
		#raw_input()
		coupling = 0
		for x in range(0,numPart):
			for y in range(x,numPart):
				if x == y:
					continue
				else:
					pass
				Rmag = math.hypot(Loc[x][0]-Loc[y][0], Loc[x][1]-Loc[y][1])
				unit_vector = (Loc[x] - Loc[y])/Rmag
				unit_dyad_term = np.dot(vec[x],vec[y])
				n_dyad_term = np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
				r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
				r_squared = (alpha*wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = (alpha*wavenumber**2)/(Rmag)
				exponent = np.exp(1j*wavenumber*Rmag)
				coupling += -(hbar/elec)*wsp*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)



		#modes[mode].append(eigen[(2*numPart)-(mode+1)]) (hbar/elec)*wsp*np.sqrt
		#print coupling
		if np.isclose(abs(np.sum(vec)),0):
			print "north north"
			modes[0].append(eigen[(2*numPart)-(mode+1)])
			interaction[0].append(coupling)
		else:
			print "north south"
			modes[1].append(eigen[(2*numPart)-(mode+1)])
			interaction[1].append(coupling)

epsb = np.linspace(1,3,21)
plt.figure()
plt.plot(epsb,modes[0],epsb,modes[1],epsb,omegas,linewidth=3)
plt.ylabel('Energy (eV)')
plt.xlabel('Embedding Medium')		
plt.legend(['NN','NS'])
#plt.show()
plt.savefig('twomer_20_withwsp.pdf')


plt.figure()
plt.plot(epsb, np.subtract(modes[0],modes[1]),epsb,np.zeros((len(epsb),1)), linewidth=3)
plt.xlabel('dielectric')
plt.ylabel('Energy (eV)')
#plt.show()
plt.savefig('NN-NS.pdf')
'''plt.figure()
plt.plot(epsb,interaction[0],epsb,interaction[1],linewidth=3)	
plt.ylabel('Energy (eV)')
plt.xlabel('Embedding Medium')	
plt.legend(['NN','NS'])
plt.show()'''
#plt.savefig('twomer_interactions.pdf')


'''New section'''

