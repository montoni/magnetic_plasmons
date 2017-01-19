import numpy as np
import math
import scipy.linalg
import scipy.sparse.linalg
import matplotlib.pyplot as plt

all_NS_vec = np.loadtxt('all_NS.txt')
out_N_in_S_vec = np.loadtxt('out_N_in_S.txt')
dipole_1_vec = np.loadtxt('dipole_1.txt')
dipole_2_vec = np.loadtxt('dipole_2.txt')
no_dipole_1_vec = np.loadtxt('no_dipole_1.txt')
no_dipole_2_vec = np.loadtxt('no_dipole_2.txt')
all_N_vec = np.loadtxt('all_N.txt')
dipole = []
all_N = []
alt_NS = []
out_N_in_S = []
no_dipole = []
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
modes = np.zeros(7,dtype=object)
eigen_unsort = np.zeros(7,dtype=object)
epsb = 1
for dist in range (90,131): # separation distance
	index = dist - 90 
	all_mag_dipole = np.zeros(7,dtype=object)
	
	center = np.zeros(7,dtype=object)
	for r in range(20,21):
		elec = 1.60217662e-19 # regular coulombs
		numPart = 24; #number of particles
		a0 = r*10**-7; #sphere radius in cm

		''' now determine geometry.'''

		# make unit vectors in centimeters.
		e1 = float(1)/float(2) * 2*(a0+(dist*.1)*10**-7) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * 2*(a0+(dist*.1)*10**-7); #long side of 30-60-90 triangle
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
		ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1 (stat-coulombs)
		hbar = 1.054571726e-34; # modified Planck in J*s
		c = 2.99e10; # speed of light in cm/s
		eps0 = 8.85418782e-12; # permittivity of free space
		epsinf = 3.77; # does this have units?
		'''Properties for silver.'''
		Eplasma = 1.46599161e-18; # J
		gamma = 0.05*elec/(hbar*16)
		wplasma = Eplasma/hbar; # plasma frequency (rad/s)
		wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		'''initialize w_0 and eigen'''
		w_0 = 3*elec/hbar
		eigen = 4*np.ones(2*numPart)
		count = 1
		for mode in range(0,10):
			while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-mode-1])) > 0.000000001:
				#print eigen[2*numPart-mode-1]
				if count == 1:
					w_0 = 0
					count = count + 1
					wsp = wsp_0
				else:
					w_0 = eigen[2*numPart-mode-1]*elec/hbar
					count = count + 1
				#print count
				alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
				msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
				tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
				gamma_ret = gamma+(tau*(wsp**2))
				gamma_eV = gamma_ret*hbar/elec
				wsp = np.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
				for n in range (0,2*numPart):
					for m in range (n,2*numPart):
						if m == n: #if m and n are the same, the hammy gets the plasmon energy
							H[n,m] = (hbar*wsp/elec)**2
						elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
							H[n,m] = 0
						else: # all other dipoles are fair game
							R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, compute the distance between them
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
				#print eigen
				print mode
			vec = np.reshape(eigenVectors[:,2*numPart-(mode+1)],(numPart,2))
			#print np.sqrt(np.sum(abs(np.dot(vec.T,all_NS_vec))))
			#print np.sqrt(np.sum(abs(np.dot(vec.T,all_N_vec))))
			#print np.sqrt(np.sum(abs(np.dot(vec.T,out_N_in_S_vec))))
			#print mode+1
			#raw_input("press enter")
			if np.sqrt(abs(np.sum(np.dot(vec.T,all_NS_vec)))) > 0.95:
				alt_NS.append(eigen[2*numPart-(mode+1)]) 
			if np.sqrt(abs(np.sum(np.dot(vec.T,out_N_in_S_vec)))) > 0.825:
				out_N_in_S.append(eigen[2*numPart-(mode+1)])
			elif np.sqrt(abs(np.sum(np.dot(vec.T,all_N_vec)))) > 0.5:
				all_N.append(eigen[2*numPart-(mode+1)])
				'''x,y = zip(*Loc)
				u,v = zip(*vec)
				plt.title(''.join(['mode = ',str(mode+1)]))
				plt.quiver(x,y,u,v)
				plt.show()'''
		if dist*.1 < 9.9:
			if len(all_N) == index+4:
				del all_N[index+3]
				del all_N[index+2]
				del all_N[index+1]
			if len(all_N) == index+3:
				del all_N[index+2]
				del all_N[index+1]
			if len(all_N) == index+2:
				del all_N[index+1]
				#del all_N[index]
		elif dist*.1 == 9.9:
			print all_N
			print len(all_N)
			print index
			raw_input("????")
		elif 13 > dist*.1 > 9.9:
			if len(all_N) == index+4:
				print all_N
				#raw_input("index+4")
				del all_N[index+3]
				del all_N[index+1]
				del all_N[index]
				del all_N[index-1]
			if len(all_N) == index+3:
				print all_N
				#raw_input("index+3")
				del all_N[index+2]
				del all_N[index]
				del all_N[index-2]
			if len(all_N) == index+2:
				print all_N
				#raw_input("index+2")
				del all_N[index+1]
				del all_N[index-1]
			if len(all_N) == index+1:
				del all_N[index-1]
		else:
			print all_N
			print len(all_N)
			print index
			raw_input("this is the end")
			#del all_N[index-1]
		if len(out_N_in_S) > index+1:
			del out_N_in_S[index+1]
		'''if len(all_N) > index+1:
			del all_N[index+1]'''
		'''elif dist*.1 == 10.5:
			if len(all_N) > index+2:
				del all_N[index+2]
				del all_N[index]'''
		#elif np.sqrt(np.sum(np.dot(vec,all_NS_vec)))
#print crossing
print len(dipole)
print len(no_dipole)
print len(alt_NS)
print len(out_N_in_S)
print len(all_N)
r = np.linspace(9,13,41)
plt.plot(r,alt_NS,r,out_N_in_S,r,all_N,linewidth=3)
plt.show()

np.savetxt('dipole_eigen',dipole)
np.savetxt('nodipole_eigen',no_dipole)
np.savetxt('alt_NS_eigen',alt_NS)
np.savetxt('all_N_eigen',all_N)
np.savetxt('out_N_in_S_eigen',out_N_in_S)
#np.savetxt('crossing',crossing)