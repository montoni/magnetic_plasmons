import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import spline


omegas = [2.76,2.755,2.745,2.74,2.735,2.73,2.715,2.71,2.7,2.685,2.67,2.655,2.64,2.625,2.595,
		  2.58,2.55,2.52,2.505,2.475,2.46,2.43,2.4,2.385,2.355,2.34,2.31,2.28,2.265,2.235]
rad = np.linspace(1,30,30)
radrad = np.linspace(1,30,291) #300 represents number of points to make between T.min and T.max

omegas_smooth = spline(rad,omegas,radrad)
NN = []
NS = []
interaction = [[],[]]
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
'''Begin by choosing a number of particles and a number of rings, defining constants and material properties'''
numPart = 12
numRings = 2
Q = [[1,0],[0,1]]
#mie_omegas = np.loadtxt('mie_omegas_eV.txt')
# [.5,-math.sqrt(3)/2],[.5,math.sqrt(3)/2],[1,0],[-.5,-math.sqrt(3)/2],[-.5,math.sqrt(3)/2],[-1,0]
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

NN_vec = np.loadtxt('hex_vec_mode_0')
NS_vec = np.loadtxt('hex_vec_mode_1')
vector = [NN_vec,NS_vec]
# create a loop over length scales, defined by particle radius
radius = np.linspace(1,30,30)
for rad in radius:
	
	a0 = rad * nm
	minor = a0
	major = 2.5*a0
	rij = 3*a0
	inter_particle_dist = 3 * a0
	theta = np.linspace(0,2*math.pi,numPart+1) # angles for placing particles around a circle
	phi = math.pi/numPart # angle for determining distance from center of circle
	dist_to_center = inter_particle_dist/(2*math.sin(phi)) # compute distance to center of ring
	ex = 3*minor/2
	ey = 3*minor*math.sqrt(3)/2
	Loc = [[-5*ex,0],[-4*ex,ey],[-4*ex,-ey],[-2*ex,ey],[-2*ex,-ey],[-ex,0],[ex,0],[2*ex,-ey],[2*ex,ey],[4*ex,-ey],[4*ex,ey],[5*ex,0]]
	'''center = 
				# quick loop to place particles in the world
				Loc = []
				for part in range(numPart):
					Loc.append(np.array([dist_to_center*math.cos(theta[part]), dist_to_center*math.sin(theta[part])]))'''
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
	alphasp = 0
	for direction in [0,1]:
		alphasp += (volume_ellipsoid/(4*math.pi))*(epsinf - epsb)/(epsinf + L_gamma[0]*(epsinf - epsb))

	omega = omegas_smooth[int((rad-1)*10)]#np.sqrt((plasma_frequency/math.sqrt(epsinf + L_gamma[0]*(epsinf - epsb)))**2 - (gamma/2)**2)
	#print omega_sp
	#raw_input()
	# initialize loop over modes, frequencies, etc.
	# this will have to be adjustable for each specific case

	count = 1
	
	w_mode = 0
	coupling = 0
	for mode in [1,2]:
		vec = vector[mode-1]
		w_mode = 0
		coupling = 0
		while abs(w_mode - (omega + np.real(coupling))) > 0.001:
			w_mode = omega + np.real(coupling)
			wavenumber = np.sqrt(epsb)*w_mode*elec/(light*hbar)
			alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
			coupling = 0
			nearfield = 0
			midfield = 0
			farfield = 0
			for x in range(0,numPart):
				for y in range(x,numPart):
					if x == y:
						continue
					else:
						pass
					R = np.subtract(Loc[(x)],Loc[(y)]) #pick out the location of each dipole, compute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					unit_vector = R/Rmag #compute unit vector between dipoles
					
					unit_dyad_term = epsb*np.dot(vec[x],vec[y])
					n_dyad_term = epsb*np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
					r_cubed = 1/(Rmag**3) #this is the 1/r^3 term (static)
					r_squared = (wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (wavenumber**2)/(Rmag)
					exponent = np.exp(1j*wavenumber*Rmag)
					coupling += -(alpha*omega)*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
					nearfield += -(alpha*omega)*r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
					midfield += (alpha*omega)*1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
					farfield += -(alpha*omega)*r_unit * exponent * (unit_dyad_term - n_dyad_term)
			
		if mode == 1:
			NN.append(w_mode)
		else:
			NS.append(w_mode)
print NN
print NS
r = np.linspace(1,30,30)
plt.figure()
plt.plot(r,NN,r,NS)
plt.legend(['NN','NS'])
plt.show()