import numpy as np
import math
import matplotlib.pyplot as plt

'''Begin by choosing a number of particles and a number of rings, defining constants and material properties'''
numPart = 3
numRings = 1
<<<<<<< HEAD
mie_omegas = np.loadtxt('mie_omegas_BEM_drude.txt')
=======
mie_omegas = np.loadtxt('mie_omegas_eV.txt')
>>>>>>> 19c0a72f65e2595fb39ed61ce4acaaf25e1ed397
c = 3.0e8 # speed of light in m/s
hbar = 1.054e-34 # hbar in J*s
nm = 1e-9 # how many nanometers are in a meter?
vacuum_perm = 8.854e-12 # Farads/meter (ugh)

elec = 1.602e-19 # electric charge (coulombs)

# properties for : silver
plasma_frequency = 9.22 # eV
gamma = 0.05 # eV, reduced by a factor of 16
epsinf = 4.98

epsb = 1
normal_modes = [[],[],[],[]]
# create a loop over length scales, defined by particle radius
radius = np.linspace(1,30,30)
for rad in radius:
<<<<<<< HEAD
	omega_sp = np.sqrt((plasma_frequency/math.sqrt(epsinf + 2*epsb))**2 - (gamma/2)**2)
	#omega_sp = mie_omegas[(rad-1)*10]
=======
	#omega_sp = plasma_frequency/math.sqrt(epsinf + 2*epsb)
	omega_sp = mie_omegas[(rad-1)*10]
>>>>>>> 19c0a72f65e2595fb39ed61ce4acaaf25e1ed397
	a0 = rad * nm
	inter_particle_dist = 2.2 * a0
	theta = np.linspace(0,2*math.pi,numPart+1) # angles for placing particles around a circle
	phi = math.pi/numPart # angle for determining distance from center of circle
	dist_to_center = inter_particle_dist/(2*math.sin(phi)) # compute distance to center of ring
	center = np.zeros((2,numRings))
	# quick loop to place particles in the world
	Loc = []
	for part in range(numPart):
		Loc.append(np.array([dist_to_center*math.cos(theta[part]), dist_to_center*math.sin(theta[part])]))
	#print Loc
	alphasp = (a0**3) * (3/(epsinf+2*epsb)) # polarizability - note: may need factor of 4*pi*eps0
	
	# initialize loop over modes, frequencies, etc.
	# this will have to be adjustable for each specific case
	H = np.zeros((2*numPart,2*numPart),dtype=float)
	eigen = np.zeros(2*numPart)
	omega_mode = np.real(omega_sp)
	count = 1
	
<<<<<<< HEAD
	for mode in range(6):
		while np.absolute(np.real(omega_mode) - np.real(eigen[2*numPart - (mode+1)])) > 0.00001:
			if count == 1:
				Q = [[1,0],[0,1],[1,0],[0,1],[1,0],[0,1]] # dipole moments in x- and y-direction
				#count = count + 1
			else:
				pass
				#Q = np.reshape(vec[2*numPart - (mode+1)],(numPart,2))
				#print Q
				#raw_input()
			omega_mode = eigen[2*numPart - (mode+1)]
=======
	for mode in range(4):
		while np.absolute(np.real(omega_mode) - np.real(eigen[2*numPart - (mode+1)])) > 0.00001:
			print omega_mode
			print eigen[2*numPart - (mode+1)]
			print omega_mode - eigen[2*numPart - (mode+1)]
			count = count +1
			print rad
			Q = [[1,0],[0,1],[1,0],[0,1],[1,0],[0,1]] # dipole moments in x- and y-direction
			omega_mode = np.real(eigen[2*numPart - (mode+1)])
>>>>>>> 19c0a72f65e2595fb39ed61ce4acaaf25e1ed397
			wavenumber = (omega_mode*elec)/(c*hbar*math.sqrt(epsb))
			alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
			mass = 1/(4*math.pi*vacuum_perm)*(elec**2)/(alpha*(omega_sp*elec/hbar)**2)
			#print mass
			#raw_input()
			for i in range(2*numPart):
				for j in range(2*numPart):
					r_ij = (Loc[i/2] - Loc[j/2])
					r_ij_mag = np.sqrt(r_ij[0]**2 + r_ij[1]**2)
					#print r_ij_mag
					nhat_ij = (Loc[i/2] - Loc[j/2])/r_ij_mag
					if i == j:
						H[i,j] = -1
					elif r_ij_mag < 0.1*nm:
						H[i,j] = 0
					else:
						#modified_distance = r_ij_mag - a0
					
						Q_dot_Q = np.dot(Q[i],Q[j])
						Q_nn_Q = np.dot(Q[i],nhat_ij)*np.dot(nhat_ij,Q[j])
							#count = count + 1
						
						# now write field terms
						near_field = (3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**3)
						int_field = -(3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**2) * 1j * wavenumber
						far_field = (Q_dot_Q - Q_nn_Q)/(r_ij_mag) * wavenumber**2
						exponential = np.exp(1j*wavenumber*r_ij_mag)
						coupling = (near_field + int_field + far_field) * exponential
						H[i,j] = -(alpha * coupling)
			zero_block = np.zeros(H.shape,dtype=complex)
			identity_block = np.identity(H.shape[0],dtype=complex)
<<<<<<< HEAD
			gamma_block = -1j*identity_block*gamma/(omega_sp)
			full_matrix = np.vstack([np.hstack([zero_block, identity_block]), np.hstack([-H, gamma_block])])
			#print full_matrix
			eigenValues, eigenVectors = np.linalg.eig(H)
			#print eigenValues
			#raw_input()
=======
			gamma_block = -identity_block*gamma
			full_matrix = np.vstack([np.hstack([zero_block, identity_block]), np.hstack([(omega_sp**2)*H, gamma_block])])
			print (omega_sp**2)*H
			eigenValues, eigenVectors = np.linalg.eig(full_matrix)
			print eigenValues
			raw_input()
>>>>>>> 19c0a72f65e2595fb39ed61ce4acaaf25e1ed397
			idx = eigenValues.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigen = np.sqrt(eigenValues[idx]) * omega_sp
			vec = eigenVectors[:,idx]
			#print vec
		#print eigen[2*numPart - (mode+1)]
		normal_modes[mode].append(eigen[2*numPart - (mode+1)])

plt.figure()
plt.plot(radius,normal_modes[0])
plt.plot(radius,normal_modes[1])
plt.plot(radius,normal_modes[2])
plt.plot(radius,normal_modes[3])
#plt.plot(radius,normal_modes[4])
#plt.plot(radius,normal_modes[5])
plt.show()




