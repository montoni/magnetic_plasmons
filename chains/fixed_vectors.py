import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline

bem_NNN = [3.614, 3.602, 3.564, 3.494, 3.412, 3.31, 3.2]
bem_N_S = [3.612, 3.601, 3.566, 3.498, 3.404, 3.294, 3.18]
bem_NSN = [3.61, 3.6, 3.57, 3.512, 3.42, 3.292, 3.12]
bem_r = [1, 5, 10, 15, 20, 25, 30]

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs
numPart = 14; #number of particles
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
mode = []
interaction = []
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
for epsb in np.linspace(0.1,5,50):
	for r in range (15,16):
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print r
		#raw_input()
		# make unit vectors in centimeters.
		rij = 3*a0
		e1 = float(1)/float(2) * (3*a0) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (3*a0); #long side of 30-60-90 triangle
		Loc = [np.array([-3*e2,2*e1]),np.array([-4*e2,e1]),np.array([-4*e2,-e1]),np.array([-3*e2,-2*e1]),
			   np.array([0, e1]),np.array([-e2, 2*e1]),np.array([-2*e2, e1]),np.array([-2*e2, -e1]),
		       np.array([-e2, -2*e1]),np.array([0, -e1]),np.array([e2, -2*e1]),np.array([2*e2, -e1]),
		       np.array([2*e2 , e1]),np.array([e2 , 2*e1])]#,np.array([3*e2,-2*e1]),np.array([4*e2,-e1]),
		       #np.array([4*e2,e1]),np.array([3*e2,2*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.
		center = [np.array([-3*e2,0]),np.array([-e2, 0]),np.array([e2,0])]#,np.array([3*e2,0])]
		'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
		H = np.zeros((2*numPart,2*numPart),dtype=complex) #initialize Hammy with zeros, twenty by twenty in this case.

		'''More constants'''
		alphasp = a0**3 * (3/(epsinf+2))
		count = 1
		#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		wsp_0 = ((mie_omegas[index])*elec/hbar)*(math.sqrt(epsinf+2))/(math.sqrt(epsinf+2*epsb))
		'''initialize w_0 and eigen'''
		w_mode = 0
		
		vec_NSN = np.loadtxt('mode_vec_0.txt')
		vec_N_S = np.loadtxt('mode_vec_1.txt')
		vec_NNN = np.loadtxt('mode_vec_2.txt')
		vec = vec_NNN
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
					nearfield += -(hbar/elec) * wsp_0/epsb * r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
					midfield += -(hbar/elec) * wsp_0/epsb * -1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
					farfield += -(hbar/elec) * wsp_0/epsb * r_unit * exponent * (unit_dyad_term - n_dyad_term)
		w_mode = wsp_0*hbar/elec+coupling
		interaction.append(coupling)
		mode.append(w_mode)
		NF[0].append(nearfield)
		IF[0].append(midfield)
		FF[0].append(farfield)

	
	#print len(NSNS)	
epsb = np.linspace(0.1,5,50)

plt.figure()
#plt.subplot(2,1,1)
#plt.plot(epsb,mode,linewidth=3)
#plt.subplot(2,1,2)
plt.plot(epsb,interaction,linewidth=3,label='NNN')
plt.scatter(epsb,np.add(NF[0],IF[0]),label = 'NNN NF + IF', color = 'C0', marker = 'o')
plt.scatter(epsb,FF[0],label = 'NNN FF', color = 'C0', marker = 's')
plt.show()
np.savetxt('dielectric_NNN.txt',np.real(mode))
np.savetxt('NNN_coupling.txt',np.real(interaction))
np.savetxt('NNN_near_mid.txt',np.real(np.add(NF[0],IF[0])))
np.savetxt('NNN_far.txt',np.real(FF[0]))


