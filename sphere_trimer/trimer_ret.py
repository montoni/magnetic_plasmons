import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm

vec =[[],[],[]]
mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
#north_north = np.loadtxt('../dielectric_study/NN_right.txt')
#north_south = np.loadtxt('../dielectric_study/NS_right.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
mag = []
<<<<<<< HEAD
for r in range(10,11):
=======
for r in np.linspace(10,10,1):
>>>>>>> 248f1ce6598f71d449290d94abddb91ea29a5b3d
	elec = 1.60217662e-19 # regular coulombs
	numPart = 3 #number of particles
	a0 = r*10**-7 #sphere radius in cm
	mie_index = (r-1)*10
	index = r

	''' now determine geometry.'''
	''' begin with electric dipole on the six-member ring '''
	# make unit vectors in centimeters.
	e1 = float(1)/float(2) * (3*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (3*a0); #long side of 30-60-90 triangle
	m_sep = 2*e2
	Loc = [np.array([0,e2]),np.array([-e1,0]),np.array([e1,0])]
		   #np.array([0,-2*e1]),np.array([e2,-e1]),np.array([e2,e1])]

	Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

	Q_mag_2 = [0,0,-1]

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
	wsp_0 = (mie_omegas[int(mie_index)])*elec/hbar
	'''initialize w_0 and eigen'''

	w_0 = wsp_0
	eigen = np.zeros(2*numPart)
<<<<<<< HEAD
	count = 1

	while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-1])) > 0.00000001:
		if count == 1:
			w_0 = 0
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
			for m in range (0,2*numPart):
				if m == n: #if m and n are the same, the hammy gets the plasmon energy
					H[n,m] = 1
					#print H[n,m]
				elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
					H[n,m] = 0
				elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
					H[m,n] = 0
				else: # all other dipoles are fair game
					R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, compute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
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
					H[n,m] = -np.real(ge) #this has the minus sign we need.
					H[m,n] = -np.real(ge)
					#print H[n,m]
		w,v = np.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
		#print w
		idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
		eigenValues = w[idx] # sorting
		eigenVectors = v[idx,:] # sorting
		eigen=np.sqrt(eigenValues)*wsp*hbar/elec # the eigenvalues have units of energy^2, so we take the square root
		print eigenVectors
		raw_input()
	print eigen[2*numPart-1]
			#print eigenVectors[0:2,2*numPart-1]
			#print eigen
		    #w_old = w_0
		    #w_0 = eigen[2*numPart-1]
		#print np.sum(eigenVectors[:,(2*numPart)-1])
		#print np.sum(eigenVectors[:,(2*numPart)-2])
		#np.savetxt('_'.join([str(r),'nm_lowest.txt']),eigenVectors[:,(2*numPart)-1])
		#np.savetxt('_'.join([str(r),'nm_second.txt']),eigenVectors[:,(2*numPart)-2])
	mag.append(eigen[2*numPart-1])


	
r = np.linspace(1,30,291)
plt.plot(r,mag,linewidth=3)	
plt.xlabel('Particle Radius (nm)')
plt.ylabel('Energy (eV)')
plt.xlim([1,30])
plt.show()
=======
>>>>>>> 248f1ce6598f71d449290d94abddb91ea29a5b3d

	count = 1
	Bfield_total = np.zeros((16,16,3),dtype=complex)
	Efield_total = np.zeros((16,16,3),dtype=complex)
	for mode in [0,1,2]:
		while np.sqrt(np.square(w_0*hbar/elec - eigen[2*numPart-(mode+1)])) > 0.00000001:
			if count == 1:
				w_0 = 0
				count = count + 1
				wsp = wsp_0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
				#if w_0 > wsp:
				#	w_0 = 3.5*elec/hbar
			#print eigen[2*numPart-1]
			wavenumber = w_0/c
			alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
			msp = (ch**2)/(alphasp*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
			#print gamma_ret
			gamma_eV = gamma_ret*hbar/elec
			wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2); # sp frequency (rad/s) corrected for radiation damping
			print wsp
			for n in range (0,2*numPart):
				for m in range (0,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1
						#print H[n,m]
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
						r_cubed = alphasp*Rmag**-3 #this is the 1/r^3 term (static)
						r_squared = (alphasp*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alphasp*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						ge = (r_unit *space_cos* (p_dot_p - p_nn_p) + (r_cubed*space_cos + r_squared*space_sin) * (3*p_nn_p - p_dot_p)) #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -(ge) #this has the minus sign we need.
						#print H[n,m]
			w,v = np.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			#print w
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting
			eigen=np.sqrt(eigenValues)*wsp*hbar/elec # the eigenvalues have units of energy^2, so we take the square root
			print eigenVectors
		vec[mode] = np.reshape(eigenVectors[:,2*numPart - (mode+1)],[numPart,2])
		
		#print vec
		#raw_input()
		'''screen_dist = 1000*a0
		
		for number in range(0,numPart):
			location = list(Loc[number])
			location.append(0)
			
			Bfield = np.empty((16,16,3),dtype=complex)
			Efield = np.empty((16,16,3),dtype=complex)
			phi_count = 0
			
			for phi in np.linspace(0,2*math.pi,16):
				theta_count = 0
				for theta in np.linspace(0,math.pi,16):
					nhat = [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]
					point = np.multiply(screen_dist,nhat)
					rmag = np.sqrt((point[0]-location[0])**2 + (point[1]-location[1])**2 + point[2]**2)
					nhat_dip = (location-point)/rmag
					Bfield[theta_count,phi_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)
					Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)
					theta_count += 1
				phi_count += 1
			Bfield_total = Bfield_total + Bfield
			Efield_total = Efield_total + Efield'''
	vec_plus = (vec[0] + vec[1])/math.sqrt(2)
	vec_minus = (vec[0] - vec[1])/math.sqrt(2)
	vec_plus_plus = (vec_plus + vec_minus)/math.sqrt(2)
	vec_minus_minus = (vec_plus - vec_minus)/math.sqrt(2)
	print vec[0]*math.sqrt(3)#_plus_plus
	print vec[1]#_minus_minus
	print vec[2]#_minus_minus
	raw_input()
	poynting = np.empty((16,16),dtype=float)
	for idx1 in range(16):
		for idx2 in range(16):
			poynting[idx1,idx2] = np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))
	theta = np.linspace(0,math.pi,16)
	phi = np.linspace(0,2*math.pi,16)
	PHI,THETA = np.meshgrid(phi,theta)
	X = poynting * np.cos(PHI)*np.sin(THETA)
	Y = poynting * np.sin(PHI)*np.sin(THETA)
	Z = poynting * np.cos(THETA)
	norm = poynting/poynting.max()
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1, projection='3d')
	plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm),linewidth=0, antialiased=False, alpha=0.5)
	plt.show()
'''np.savetxt('collinear',coll)
np.savetxt('anti_collinear',anti_coll)
np.savetxt('parallel',para)
np.savetxt('anti_parallel',anti_para)'''