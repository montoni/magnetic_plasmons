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
mode = []
interaction = []
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
for r in range(1,31):
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
	alphasp = a0**3 * (3/(epsinf+2))
	count = 1
	#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	wsp_0 = ((mie_omegas[index])*elec/hbar)#*(math.sqrt(epsinf+2))/(math.sqrt(epsinf+2*epsb))
	'''initialize w_0 and eigen'''
	w_mode = 0
	
	Alt = np.loadtxt('coro_0_vec')
	RadNode = np.loadtxt('coro_1_vec')
	TwoNode = np.loadtxt('coro_2_vec')
	Dipole = np.loadtxt('coro_4_vec')
	All_N = np.loadtxt('coro_6_vec')
	#Alt = np.loadtxt('coro_0_vec.txt')
	#Alt = np.loadtxt('coro_0_vec.txt')
	vec = All_N
	coupling = 0
	count = 0
	while abs(np.real(w_mode)*hbar/elec - (np.real(wsp_0)*hbar/elec + np.real(coupling))) > 0.000001:
		if count == 0:
			w_mode = 0
			count += 1
			print count
		else:
			w_mode = np.real((wsp_0*hbar/elec + coupling) * elec/hbar)
			count += 1
			#print count
			print w_mode*hbar/elec
		wavenumber = math.sqrt(epsb)*w_mode/c
		alpha = ((alphasp**-1) - 1j*(2./3.)*wavenumber**3)**-1
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
				Rmag = math.hypot(Loc[x][0]-Loc[y][0], Loc[x][1]-Loc[y][1])
				unit_vector = (Loc[x] - Loc[y])/Rmag
				unit_dyad_term = np.dot(vec[x],vec[y])
				n_dyad_term = np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
				r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
				r_squared = (alpha*wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = (alpha*wavenumber**2)/(Rmag)
				exponent = np.exp(1j*wavenumber*Rmag)
				coupling += -(hbar/elec)*wsp_0/epsb*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
				nearfield += -(hbar/elec) * wsp_0 * r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
				midfield += -(hbar/elec) * wsp_0 * -1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
				farfield += -(hbar/elec) * wsp_0 * r_unit * exponent * (unit_dyad_term - n_dyad_term)
	w_mode = np.real(wsp_0*hbar/elec+coupling)
	interaction.append(coupling)
	mode.append(w_mode)
	NF[0].append(nearfield)
	IF[0].append(midfield)
	FF[0].append(farfield)


#print len(NSNS)	
epsb = np.linspace(1,30,30)

plt.figure()
plt.subplot(2,1,1)
plt.plot(epsb,mode,linewidth=3)
plt.subplot(2,1,2)
plt.plot(epsb,interaction,linewidth=3)
plt.scatter(epsb,np.add(NF[0],IF[0]),label = 'NF + IF', color = 'C0', marker = 'o')
plt.scatter(epsb,FF[0],label = 'FF', color = 'C0', marker = 's')
plt.show()
#np.savetxt('dielectric_N_S.txt',np.real(mode))
np.savetxt('allN_coupling.txt',np.real(interaction))
np.savetxt('allN_near_mid.txt',np.real(np.add(NF[0],IF[0])))
np.savetxt('allN_far.txt',np.real(FF[0]))


