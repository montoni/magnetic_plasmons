import numpy as np
import math
import matplotlib.pyplot as plt

'''Begin by choosing a number of particles and a number of rings, defining constants and material properties'''
numPart = 2
numRings = 1
mie_omegas = np.loadtxt('mie_omegas_eV.txt')
c = 3.0e8 # speed of light in m/s
hbar = 1.054e-34 # hbar in J*s
nm = 1e-9 # how many nanometers are in a meter?
vacuum_perm = 8.854e-12 # Farads/meter (ugh)

elec = 1.602e-19 # electric charge (coulombs)

# properties for : silver
plasma_frequency = 9.15 # eV
gamma = 0.05/16 # eV, reduced by a factor of 16
epsinf = 3.77

epsb = 1
normal_modes = [[],[],[],[]]
# create a loop over length scales, defined by particle radius
radius = np.linspace(1,30,30)
for rad in radius:
	#omega_sp = plasma_frequency/math.sqrt(epsinf + 2*epsb)
	omega_sp = mie_omegas[(rad-1)*10]
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
	print Loc
	alphasp = (a0**3) * (3/(epsinf+2*epsb)) # polarizability - note: may need factor of 4*pi*eps0
	
	# initialize loop over modes, frequencies, etc.
	# this will have to be adjustable for each specific case
	H = np.zeros((2*numPart,2*numPart),dtype=complex)
	eigen = np.zeros(2*numPart)
	omega_mode = omega_sp
	count = 1
	
	for mode in range(4):
		while np.absolute(omega_mode - eigen[2*numPart - (mode+1)]) > 0.00001:
			if count == 1:
				Q = [[1,0],[0,1],[1,0],[0,1]] # dipole moments in x- and y-direction
				#count = count + 1
			else:
				pass
				#Q = np.reshape(vec[2*numPart - (mode+1)],(numPart,2))
				#print Q
				#raw_input()
			omega_mode = eigen[2*numPart - (mode+1)]
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
						H[i,j] = 1
					elif r_ij_mag < 0.1*nm:
						H[i,j] = 0
					else:
						#modified_distance = r_ij_mag - a0
						if count == 1:
							Q_dot_Q = np.dot(Q[i],Q[j])
							Q_nn_Q = np.dot(Q[i],nhat_ij)*np.dot(nhat_ij,Q[j])
							#count = count + 1
						else:
							Q_dot_Q = np.dot(Q[i/2],Q[j/2])
							Q_nn_Q = np.dot(Q[i/2],nhat_ij)*np.dot(nhat_ij,Q[j/2])
						# now write field terms
						near_field = (3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**3)
						int_field = -(3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**2) * 1j * wavenumber
						far_field = 2*(Q_dot_Q - Q_nn_Q)/(r_ij_mag) * wavenumber**2
						exponential = np.exp(1j*wavenumber*r_ij_mag)
						coupling = (near_field + int_field + far_field) * exponential
						H[i,j] = -(alpha * coupling)
			zero_block = np.zeros(H.shape,dtype=complex)
			identity_block = np.identity(H.shape[0],dtype=complex)
			gamma_block = -1j*identity_block*gamma/(omega_sp)
			full_matrix = np.vstack([np.hstack([zero_block, identity_block]), np.hstack([-H, gamma_block])])
			#print full_matrix
			eigenValues, eigenVectors = np.linalg.eig(H)
			#print eigenValues
			#raw_input()
			idx = eigenValues.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigen = np.sqrt(eigenValues[idx]) * omega_sp
			vec = eigenVectors[:,idx]
			#print vec
		print eigen[2*numPart - (mode+1)]
		normal_modes[mode].append(eigen[2*numPart - (mode+1)])

plt.figure()
plt.plot(radius,normal_modes[0])
plt.plot(radius,normal_modes[1])
plt.plot(radius,normal_modes[2])
plt.plot(radius,normal_modes[3])
#plt.plot(radius,normal_modes[4])
#plt.plot(radius,normal_modes[5])
plt.show()




