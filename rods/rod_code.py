import numpy as np
import math
import matplotlib.pyplot as plt

'''Begin by choosing a number of particles and a number of rings, defining constants and material properties'''
numPart = 3
numRings = 1

#mie_omegas = np.loadtxt('mie_omegas_eV.txt')

light = 3.0e8 # speed of light in m/s
hbar = 1.054e-34 # hbar in J*s
nm = 1e-9 # how many nanometers are in a meter?
vacuum_perm = 8.854e-12 # Farads/meter (ugh)

elec = 1.602e-19 # electric charge (coulombs)

# properties for : silver
plasma_frequency = 9.15 # eV
gamma = 0.05 # eV, reduced by a factor of 16
epsinf = 3.77

epsb = 1
normal_modes = [[],[],[],[],[],[]]
# create a loop over length scales, defined by particle radius
radius = np.linspace(1,40,40)
for rad in radius:
	
	a0 = rad * nm
	inter_particle_dist = 3 * a0
	theta = np.linspace(0,2*math.pi,numPart+1) # angles for placing particles around a circle
	phi = math.pi/numPart # angle for determining distance from center of circle
	dist_to_center = inter_particle_dist/(2*math.sin(phi)) # compute distance to center of ring
	center = np.zeros((2,numRings))
	# quick loop to place particles in the world
	Loc = []
	for part in range(numPart):
		Loc.append(np.array([dist_to_center*math.cos(theta[part]), dist_to_center*math.sin(theta[part])]))
	#print Loc
	#raw_input()
	# particle axes lengths
	a = 2.5*a0
	b = a0
	c = a0
	volume_ellipsoid = (4./3.)*math.pi*a*b*c
	if a > b:
		pro_e = math.sqrt(1 - (b/a)**2)
		L_x = ((1 - pro_e**2)/(2*pro_e**3))*(math.log((1 + pro_e)/(1 - pro_e)) - 2*pro_e)
		L_yz = 0.5*(1 - L_x)
		L_gamma = [L_x,L_yz]
	# oblate
	if a < b:
		obl_e = math.sqrt((a/c)**2 - 1)
		L_x = ((1 + obl_e**2)/(obl_e**3))*(obl_e - math.arctan(obl_e))
		L_yz = 0.5*(1 - L_z)
		L_gamma = [L_x,L_yz]
	if a == b:
		L_gamma = [1/3,1/3]
	# generate polarizability

	for direction in [0,1]:
		alphasp = (volume_ellipsoid/(4*math.pi))*(epsinf - epsb)/(epsinf + L_gamma[0]*(epsinf - epsb))

	#print alphasp
	#raw_input()
	omega_sp = 2.8#np.sqrt((plasma_frequency/math.sqrt(epsinf + L_gamma[0]*(epsinf - epsb)))**2 - (gamma/2)**2)
	#print omega_sp
	#raw_input()
	# initialize loop over modes, frequencies, etc.
	# this will have to be adjustable for each specific case
	H = np.zeros((2*numPart,2*numPart),dtype=float)
	eigen = np.zeros(2*numPart)
	omega_mode = omega_sp
	count = 1

	for mode in range(2*numPart):
		while np.absolute((omega_mode) - (eigen[2*numPart - (mode+1)])) > 0.0001:
			Q = [[1,0],[0,1]] # dipole moments in x- and y-direction
			
			omega_mode = (eigen[2*numPart - (mode+1)])
			print omega_mode

			wavenumber = (omega_mode*elec)/(light*hbar*math.sqrt(epsb))
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
						Q_dot_Q = np.dot(Q[i%2],Q[j%2])
						Q_nn_Q = np.dot(Q[i%2],nhat_ij)*np.dot(nhat_ij,Q[j%2])
						# now write field terms
						near_field = (3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**3)
						int_field = -(3*Q_nn_Q - Q_dot_Q)/(r_ij_mag**2) * 1j * wavenumber
						far_field = (Q_dot_Q - Q_nn_Q)/(r_ij_mag) * wavenumber**2
						exponential = np.exp(1j*wavenumber*r_ij_mag)
						coupling = (near_field + int_field + far_field) * exponential
						H[i,j] = -(alpha * coupling)
			#print H
			#raw_input()
			v, w = np.linalg.eig(H)
			idx = v.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigen = np.sqrt(v[idx]) * omega_sp
			eigenVectors = w[:,idx]
		normal_modes[mode].append(eigen[2*numPart - (mode+1)])
		vec = np.reshape(eigenVectors[:,2*numPart - (mode+1)],(numPart,2))
		x,y = zip(*Loc)
		u,v = zip(*vec)
		'''plt.figure()
		plt.quiver(x,y,u,v)
		plt.show()'''


plt.figure()
plt.plot(radius,normal_modes[0])
plt.plot(radius,normal_modes[1])
plt.plot(radius,normal_modes[2])
plt.show()




