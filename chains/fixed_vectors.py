import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm

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
gamma = 0.05*elec/(hbar)
NNN = []
NSN = []
N_S = []
interaction = [[],[],[]]
NF = [[],[],[]]
IF = [[],[],[]]
FF = [[],[],[]]
epsb = 1.77
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
vec_1 = np.loadtxt('mode_vec_0.txt')
vec_2 = np.loadtxt('mode_vec_1.txt')
vec_3 = np.loadtxt('mode_vec_2.txt')
vectors = [[vec_1],[vec_2],[vec_3]]
for sep in np.linspace(1,1,1):
	for r in range (1,31):
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print r
		#raw_input()
		# make unit vectors in centimeters.
		rij = (2+sep)*a0
		e1 = float(1)/float(2) * (rij) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (rij); #long side of 30-60-90 triangle
		Loc = [np.array([-3*e2,2*e1]),np.array([-4*e2,e1]),np.array([-4*e2,-e1]),np.array([-3*e2,-2*e1]),
			   np.array([0, e1]),np.array([-e2, 2*e1]),np.array([-2*e2, e1]),np.array([-2*e2, -e1]),
		       np.array([-e2, -2*e1]),np.array([0, -e1]),np.array([e2, -2*e1]),np.array([2*e2, -e1]),
		       np.array([2*e2 , e1]),np.array([e2 , 2*e1])]
		center = [np.array([-3*e2,0]),np.array([-e2, 0]),np.array([e2,0])]
		'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
		H = np.zeros((2*numPart,2*numPart),dtype=complex) #initialize Hammy with zeros, twenty by twenty in this case.

		'''More constants'''
		alphasp = a0**3 * (3/(epsinf+2*epsb))
		count = 1
		#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		wsp_0 = ((mie_omegas[index]-r*0.0034)*elec/hbar)*(math.sqrt(epsinf+2))/(math.sqrt(epsinf+2*epsb))
		'''initialize w_0 and eigen'''
		w_mode = 0
		
		

		coupling = 0
		count = 0
		#Bfield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		#Efield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		for mode in [0,1,2]:
			vec = np.reshape(vectors[mode],(numPart,2))

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
				alpha = ((alphasp**-1) - 1j*(2./(3.*epsb))*wavenumber**3)**-1
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
						unit_dyad_term = epsb*np.dot(vec[x],vec[y])
						n_dyad_term = epsb*np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*wavenumber**2)/(Rmag)
						exponent = np.exp(1j*wavenumber*Rmag)
						coupling += -(hbar/elec)*wsp_0*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
						nearfield += -(hbar/elec) * wsp_0 * r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
						midfield += -(hbar/elec) * wsp_0 * -1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
						farfield += -(hbar/elec) * wsp_0 * r_unit * exponent * (unit_dyad_term - n_dyad_term)
			w_mode = wsp_0*hbar/elec + np.real(coupling)
			if mode == 0:
				NSN.append(w_mode)
				NF[1].append(nearfield)
				IF[1].append(midfield)
				FF[1].append(farfield)
				interaction[1].append(coupling)
			elif mode == 1:
				N_S.append(w_mode)
				NF[2].append(nearfield)
				IF[2].append(midfield)
				FF[2].append(farfield)
				interaction[2].append(coupling)
			else:
				NNN.append(w_mode)
				NF[0].append(nearfield)
				IF[0].append(midfield)
				FF[0].append(farfield)
				interaction[0].append(coupling)
r = np.linspace(1,30,30)
'''np.savetxt('water_NNN_30_space.txt',NNN)
np.savetxt('water_NSN_30_space.txt',NSN)
np.savetxt('water_N_S_30_space.txt',N_S)'''
plt.figure()
plt.plot(r,NNN,r,NSN,r,N_S)
plt.savefig('water_linear_scale_eig.pdf')
plt.show()

plt.figure()
plt.plot(r,interaction[0],label='NNN',linewidth=3)
plt.plot(r,interaction[1],label='NSN',linewidth=3)
plt.plot(r,interaction[2],label='N_S',linewidth=3)
plt.scatter(r,np.add(NF[0],IF[0]),label='NNN near+mid',color='C0',marker='o')
plt.scatter(r,np.add(NF[1],IF[1]),label='NSN near+mid',color='C1',marker='o')
plt.scatter(r,np.add(NF[2],IF[2]),label='N_S near+mid',color='C2',marker='o')
plt.scatter(r,FF[0],label='NNN far',color='C0',marker='s')
plt.scatter(r,FF[1],label='NSN far',color='C1',marker='s')
plt.scatter(r,FF[2],label='N_S far',color='C2',marker='s')
plt.ylabel('Energy (eV)')
plt.xlabel('Particle Radius r_0 (nm)')	
plt.legend()
plt.savefig('water_linear_scale_int.pdf')
plt.show()
