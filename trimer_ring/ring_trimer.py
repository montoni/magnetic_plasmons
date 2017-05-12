import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline

bem_r = [1, 5, 10, 15, 20, 25]
bem_NN = [3.616, 3.604, 3.558, 3.482, 3.402, 3.31]
bem_NS = [3.61, 3.6, 3.568, 3.51, 3.416, 3.288]

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
interaction = [[],[]]
elec = 1.60217662e-19 # regular coulombs
numPart = 13; #number of particles
me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; # does this have units?
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar*16)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
for epsb in range(1,2):
	NN = []
	NS = []
	for r in range (1,31):
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print r
		# make unit vectors in centimeters.
		rij = 3*a0
		e1 = float(1)/float(2) * (3*a0) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (3*a0); #long side of 30-60-90 triangle
		Loc = [np.array([0, e1]),np.array([-e2, 2*e1]),
			   np.array([-2*e2, e1]),np.array([-2*e2, -e1]),
			   np.array([-e2, -2*e1]),np.array([0, -e1]),
			   np.array([e2, -2*e1]),np.array([2*e2, -e1]),
			   np.array([2*e2 , e1]),np.array([e2 , 2*e1]),
			   np.array([-e2,4*e1]),np.array([0,5*e1]),np.array([e2,4*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.

		'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
		H = np.zeros((2*numPart,2*numPart),dtype=complex) #initialize Hammy with zeros, twenty by twenty in this case.

		'''More constants'''

		count = 1
		#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		wsp_0 = (mie_omegas[index])*elec/hbar
		'''initialize w_0 and eigen'''
		w_0 = 3*elec/hbar
		eigen = np.ones(2*numPart)
		for mode in range(0,3):
			while np.sqrt(np.square(w_0*hbar/elec - eigen[(2*numPart)-(mode+1)])) > 0.00000001:
				if count == 1:
					wsp = wsp_0
					count = count + 1
					w_0 = 0
				else:
					count = count + 1
					w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
				alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
				wavenumber = (w_0)/(c*math.sqrt(epsb))
				alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
				msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
				tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
				gamma_ret = gamma+tau*(wsp**2) # I pulled this from the beats paper
				gamma_eV = gamma_ret*hbar/elec
				#wsp = wsp_0
				wsp = wsp_0 # sp frequency (rad/s) corrected for radiation damping
				for n in range (0,2*numPart):
					for m in range (0,2*numPart):
						if m == n: #if m and n are the same, the hammy gets the plasmon energy
							H[n,m] = 1
						elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
							H[n,m] = 0
						elif n == m+1 and m%2 == 0:
							H[m,n] = 0
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
							H[n,m] = -ge #this has the minus sign we need.
				w,v = np.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
				idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
				eigenValues = w[idx] # sorting
				eigenVectors = v[:,idx] # sorting
				eigen=np.sqrt(eigenValues) * wsp * hbar/elec # the eigenvalues have units of energy^2, so we take the square root
				#print eigen
			    #w_old = w_0
			    #w_0 = eigen[2*numPart-1]
			vec = eigenVectors[:,(2*numPart)-(mode+1)]
			vec = np.reshape(vec,[numPart,2])
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
					r_squared = (alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (alpha*w_0**2)/(Rmag*(c**2))
					exponent = np.exp(1j*w_0*Rmag/c)
					coupling += -(hbar/elec)*wsp*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
			'''x,y = zip(*Loc)
			u,v = zip(*vec)
			plt.title(''.join(['radius = ',str(r)]))
			plt.quiver(x,y,u,v)
			plt.show()
			raw_input()'''
			if abs(np.sum(vec)) < 10**-5:
				NN.append(eigen[2*numPart-mode-1])
				interaction[0].append(coupling)
			else:
				NS.append(eigen[2*numPart-mode-1])
				interaction[1].append(coupling)
	

	smooth_r = np.linspace(1,25,25)
	NN_smooth = spline(bem_r,bem_NN,smooth_r)
	NS_smooth = spline(bem_r,bem_NS,smooth_r)

	print len(NN)
	print len(NS)
	NS = np.reshape(NS,[30,2])
	NS_int = np.reshape(interaction[1],(30,2))
	r = np.linspace(1,30,30)
	plt.figure()
	plt.plot(r,NN,r,NS[:,0],smooth_r,NN_smooth,smooth_r,NS_smooth,linewidth=3)	
	plt.legend(['NN','NS','BEM NN','BEM NS'])
	plt.ylabel('Energy (eV)')
	plt.xlabel('Particle Radius r_0 (nm)')
	#plt.savefig('threemer_eigenvalues.pdf')
	plt.show()

	plt.figure()
	plt.plot(r,interaction[0],r,NS_int[:,0],linewidth=3)
	plt.legend(['NN','NS'])
	plt.ylabel('Energy (eV)')
	plt.xlabel('Particle Radius r_0 (nm)')
	#plt.savefig('threemer_interactions.pdf')
	plt.show()
	#np.savetxt('_'.join([str(epsb),'NN.txt']),NN)
	#np.savetxt('_'.join([str(epsb),'NS.txt']),NS)
#np.savetxt('crossing',crossing)
#print static_NN
#print static_NS
