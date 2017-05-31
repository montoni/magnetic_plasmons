import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline

static_NN = []
static_NS = []
bem_NN = [3.614, 3.603, 3.562, 3.5, 3.405, 3.295, 3.17]
bem_NS = [3.61, 3.6, 3.569, 3.509, 3.413, 3.276, 3.112]
mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs

numPart = 10; #number of particles

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
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
NFIF = [[],[]]
for r in np.linspace(1,1,1):
	a0 = r*10**-7; #sphere radius in cm
	alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
	index = r
	''' now determine geometry.'''
	print index
	# make unit vectors in centimeters.
	rij = 3*a0
	e1 = .5*rij
	e2 = rij*math.sqrt(3)/2
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
	
	Loc = places
	'''Loc = np.zeros((numPart,2))
	for x in places:
		print x
		if np.isclose(x.all(),Loc[:].any()):
			pass
		else:
			Loc[count] = x'''


	#Loc = np.unique(Loc)
	#print places
	
	#raw_input()
	#plt.show()
	#raw_input()
	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = (np.zeros((2*numPart,2*numPart),dtype=float))#initialize Hammy with zeros, twenty by twenty in this case.

	'''More constants'''

	count = 1
	#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	wsp_0 = (mie_omegas[index*10 - 10])*elec/hbar
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
			wavenumber = (w_0)/(c*math.sqrt(epsb))
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
						r_squared = (alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						exponent = np.exp(1j*w_0*Rmag/c)
						ge = ((r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))) * exponent #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -(ge)#*(hbar/elec)*wsp #this has the minus sign we need.
						H[m,n] = np.conj(-ge)
						

					

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
		grid_x = 2.5*e2
		grid_y = 2.5*e1
		size = 251
		Bfield_total = np.zeros((size,size),dtype=float)
		for number in range(0,10):
			point = Loc[number]
			xcount = 0
			Bfield = np.empty((size,size),dtype=float)
			for xx in np.linspace(-grid_y,grid_y,size):
				ycount = 0
				for yy in np.linspace(-grid_x,grid_x,size):
					rmag = np.sqrt((point[0]-yy)**2 + (point[1]-xx)**2)
					nhat = -([yy,xx]-point)/rmag
					#print rmag
					#raw_input()
					#vec[number][0],vec[number][1] = vec[number][1],vec[number][0]
					if rmag < .5*a0:
						Bfield[xcount,ycount] = 0
					else:
						Bfield[xcount,ycount] = np.cross(nhat,vec[number])*(rmag**-1)*(1-((1j*wavenumber*rmag)**-1))
					ycount += 1
				xcount += 1
			Bfield_total = Bfield_total + Bfield
		# ,levels=np.linspace(np.amin(Bfield_total),np.amax(Bfield_total),30)
		xx = np.linspace(-grid_x,grid_x,size)
		yy = np.linspace(-grid_y,grid_y,size)
		Loc = places
		u,v = zip(*vec)
		xloc, yloc = zip(*Loc)
		plt.figure()
		plt.contourf(xx,yy,Bfield_total,cmap='bwr',levels=np.linspace(np.amin(Bfield_total),np.amax(Bfield_total),30))
		#plt.colorbar()
		plt.tight_layout()
		plt.scatter(xloc,yloc)
		plt.quiver(xloc,yloc,u,v)
		plt.savefig('mode_'+str(mode)+ '_field.pdf')
		plt.show()