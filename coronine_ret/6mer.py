import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt


mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
dipole = [] # doubly degenerate
no_dipole = [] # doubly degenerate
all_N = [] #all north
out_N_in_S = [] # that weird "excited state"
alt_NS = [] # maximally out of phase magnets
center = np.zeros(7,dtype=object)
for r in range(1,2):
	elec = 1.60217662e-19 # regular coulombs
	numPart = 24; #number of particles
	a0 = r*10**-7; #sphere radius in cm
	index = (r - 1) * 10
	''' now determine geometry.'''

	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (3*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (3*a0); #long side of 30-60-90 triangle
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
	wsp_0 = mie_omegas[index]*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = 3*elec/hbar
	eigen = np.ones(2*numPart)
	count = 1
	for mode in range(0,7):
		while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-(mode+1)])) > 0.00000001:
			if count == 1:
				w_0 = 0
				wsp = wsp_0
				count = count + 1
			else:
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
				count = count + 1
			#print count
			#print eigen[2*numPart-(mode+1)]
			alphasp = (a0**3)*(3/(epsinf+2)); # polarizability (cm^3)
			wavenumber = (w_0)/(c*math.sqrt(epsb))
			alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
			msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(wsp**2) # I pulled this from the beats paper
			gamma_eV = gamma_ret*hbar/elec
			wsp = wsp_0#math.sqrt((wsp_0)**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,2*numPart):
				for m in range (0,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1#(hbar*wsp/elec)
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alpha/Rmag**3 #this is the 1/r^3 term (static)
						r_squared = (alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -np.real(ge) #this has the minus sign we need.
			#diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
			#Ht = np.matrix.transpose(H) # this is the transpose of H
			#Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			#Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			w,v = scipy.linalg.eigh(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting
			eigen=np.sqrt(eigenValues)*wsp*hbar/elec # the eigenvalues have units of energy^2, so we take the square root
			#print eigen
		vec = np.reshape(eigenVectors[:,2*numPart-(mode+1)],(numPart,2))
		x,y = zip(*Loc)
		u,v = zip(*vec)
		plt.title(''.join(['mode = ',str(mode+1)]))
		plt.quiver(x,y,u,v)
		plt.show()
		#np.savetxt('coro_' + str(mode) + '_vec',vec)
		'''mag_dipole = []
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
		if abs(np.sum(vec)) < 10**-10:
			if abs(np.sum(mag_dipole)) < 10**-12:
				if abs(mag_dipole[0]) < 10**-12:
					if np.dot(vec[12],vec[21]) > 0: 
						alt_NS.append(np.real(eigen[2*numPart-(mode+1)]))
					else:
						no_dipole.append(np.real(eigen[2*numPart-(mode+1)]))
				else:
					alt_NS.append
			else:
				if np.dot(vec[0],vec[21]) * np.dot(vec[0],vec[3]) < 0 :
					all_N.append(np.real(eigen[2*numPart-(mode+1)]))
				elif abs(np.sum(vec[21])) < 10**-5:
					all_N.append(np.real(eigen[2*numPart-(mode+1)]))
				elif np.dot(vec[0],vec[21]) > 0:
					all_N.append(np.real(eigen[2*numPart-(mode+1)]))
				elif mag_dipole[0] * mag_dipole[1] > 0:
					all_N.append(np.real(eigen[2*numPart-(mode+1)]))
				else:
					out_N_in_S.append(np.real(eigen[2*numPart-(mode+1)]))
		else:
			dipole.append(np.real(eigen[2*numPart-(mode+1)]))
	if len(out_N_in_S) > len(all_N):
		all_N.append(out_N_in_S[index])
		del out_N_in_S[index]
	if len(alt_NS) > len(all_N):
		all_N.append(alt_NS[index])
		del alt_NS[index]
	if len(all_N) > len(out_N_in_S):
		out_N_in_S.append(all_N[index])
		del all_N[index]
	if len(out_N_in_S) > alt_NS:
		alt_NS.append(out_N_in_S[index])
		del out_N_in_S[index]
	if len(out_N_in_S) > index + 1:
		if len(out_N_in_S) > len(no_dipole)*0.5:
			print len(out_N_in_S)
			print len(no_dipole)
			print index
			del out_N_in_S[index]
	if len(dipole)*.5 > len(alt_NS):
		alt_NS.append(dipole[index*2])
		del dipole[index*2]
	if out_N_in_S[index-1] - out_N_in_S[index] > 0.01:
		all_N.append(out_N_in_S[index])
		out_N_in_S.append(all_N[index])
		del all_N[index]
		del out_N_in_S[index]'''
	'''if index == 290:
		all_N.append(out_N_in_S[index])
		alt_NS.append(all_N[index])
		out_N_in_S.append(alt_NS[index])
		del all_N[index]
		del alt_NS[index]
		del out_N_in_S[index]'''
	'''if out_N_in_S[index-1] < all_N[index-1]:
		print index
		print len(out_N_in_S)
		print len(all_N)
		if out_N_in_S[index] > all_N[index]:
			crossing = ((index+index-1)*0.5*0.1)+1
	if len(out_N_in_S) > index+1:
		del out_N_in_S[index]
	print len(dipole)
	print len(no_dipole)
	print len(alt_NS)
	print len(out_N_in_S)
	print len(all_N)'''
	#raw_input("press enter")


'''print crossing
print len(dipole)
print len(no_dipole)
print len(alt_NS)
print len(out_N_in_S)
print len(all_N)
dipole = np.reshape(dipole,[291,2])
no_dipole = np.reshape(no_dipole,[291,2])
r = np.linspace(1,30,291)
plt.plot(r,alt_NS,r,out_N_in_S,r,all_N,r,dipole[:,0],r,no_dipole[:,0],linewidth=3)
plt.ylabel('Energy (eV)')
plt.xlabel('Particle Radius (nm) / Scale Factor')
plt.legend(['Alternating North-South','South-Inside-North-Outside','All-North','Net Electric Dipole','No Net Dipoles'],loc=3)
plt.show()

np.savetxt('dipole_eigen',dipole)
np.savetxt('nodipole_eigen',no_dipole)
np.savetxt('alt_NS_eigen',alt_NS)
np.savetxt('altl_N_eigen',all_N)
np.savetxt('out_N_in_S_eigen',out_N_in_S)
#np.savetxt('crossing',crossing)'''