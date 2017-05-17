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
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
for epsb in range(1,2):
	NNN = []
	NSN = []
	N_S = []
	interaction = [[],[],[]]
	#NSNS = []
	for r in range (1,31):
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
		wsp_0 = (mie_omegas[index])*elec/hbar
		'''initialize w_0 and eigen'''
		w_mode = wsp_0
		
		vec = [[-1,0],[-.5,-math.sqrt(3)/2],[.5,-math.sqrt(3)/2],[1,0],[-1,0],[-1,0],[-1,0],[1,0],[1,0],[1,0],[1,0],[.5,math.sqrt(3/2)],[-.5,math.sqrt(3/2)],[-1,0]]
		#vec = [[-1,0],[-.5,-math.sqrt(3)/2],[.5,-math.sqrt(3)/2],[1,0],[0,-1],[1,0],[0,1],[0,1],[-1,0],[0,-1],[1,0],[.5,math.sqrt(3/2)],[-.5,math.sqrt(3/2)],[-1,0]]
		#vec = [[-1,0],[-.5,-math.sqrt(3)/2],[.5,-math.sqrt(3)/2],[1,0],[.5,math.sqrt(3)/2],[0,0],[-.5,math.sqrt(3)/2],[.5,math.sqrt(3)/2],[0,0],[-.5,math.sqrt(3)/2],[-1,0],[-.5,-math.sqrt(3/2)],[.5,-math.sqrt(3/2)],[1,0]]
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
				unit_dyad_term = np.dot(vec[x]/math.sqrt(14),vec[y]/math.sqrt(14))
				n_dyad_term = np.dot(vec[x]/math.sqrt(14),unit_vector)*np.dot(unit_vector,vec[y]/math.sqrt(14))
				r_cubed = alphasp/(Rmag**3) #this is the 1/r^3 term (static)
				r_squared = (alphasp*w_mode)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
				r_unit = (alphasp*w_mode**2)/(Rmag*(c**2))
				exponent = np.exp(1j*w_mode*Rmag/c)
				coupling += -(hbar/elec)*wsp_0*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
		w_mode = wsp_0*hbar/elec+coupling
		NNN.append(w_mode)

	
	#print len(NSNS)	
	r = np.linspace(1,30,30)

	plt.figure()
	plt.plot(r,NNN,linewidth=3)
	plt.show()

	np.savetxt('NNN.txt',NNN)