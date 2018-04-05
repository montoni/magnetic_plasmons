import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline

static_NN = []
static_NS = []
bem_NN = [3.614, 3.603, 3.562, 3.5, 3.405, 3.295, 3.17]
bem_NS = [3.61, 3.6, 3.569, 3.509, 3.413, 3.276, 3.112]
mie_omegas = np.loadtxt('../mie/mie_omegas_vacuum.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs

numPart = 31 #number of particles

me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; 
'''Properties for silver.'''
Eplasma = 9.15*elec/hbar; # J
gamma = 0.05*elec/(hbar)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
epsb = 1
NN = []
NS = []
interaction = [[],[]]
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
NN_vec = np.loadtxt('31_0_vec.txt')
NS_vec = np.loadtxt('31_1_vec.txt')
vectors = [NN_vec,NS_vec]
for r in range(1,51):
	a0 = r*10**-7; #sphere radius in cm
	alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
	index = r-1
	''' now determine geometry.'''
	print index
	# make unit vectors in centimeters.
	rij = 3*a0#2*a0+(2.5*10**-7)
	part_per_ring = numPart
	theta = 2*math.pi/(6)
	phi = theta/2.
	#print theta
	#print phi
	#print np.cos(phi)
	#print np.sin(phi)
	cent_dist = rij/(2*np.tan(phi))
	part_to_cent = math.sqrt((rij/2)**2 + (cent_dist)**2)
	second_layer = math.sqrt(3) * rij
	third_layer = 2*rij
	fourth_layer = math.sqrt(7)*rij
	center = np.array([0,0])
	#print centers
	places = []
	places.append(center)
	if part_per_ring == 1:
		continue
	elif part_per_ring >= 7:
		for num in range(part_per_ring-1):
			if num < 6:
				places.append(center + np.array([(part_to_cent*np.cos(theta*(num))),(part_to_cent*np.sin(theta*(num)))]))
			elif 6 <= num < 12:
				places.append(center + np.array([(second_layer*np.cos(phi+theta*(num))),(second_layer*np.sin(phi+theta*(num)))]))
			elif 12 <= num < 18:
				places.append(center + np.array([(third_layer*np.cos(theta*(num))),(third_layer*np.sin(theta*(num)))]))
			elif 18 <= num < 24:
				places.append(center + np.array([(fourth_layer*np.cos(2*phi/3+theta*(num))),(fourth_layer*np.sin(2*phi/3+theta*(num)))]))
			elif 24 <= num < 30:
				places.append(center + np.array([(fourth_layer*np.cos(4*phi/3+theta*(num))),(fourth_layer*np.sin(4*phi/3+theta*(num)))]))
	
	count = 0 
	#print places
	#raw_input()


	#Loc = np.unique(Loc)
	#print places
	Loc = places
	xloc, yloc = zip(*Loc)
	#plt.scatter(xloc,yloc)
	#raw_input()
	#plt.show()
	#raw_input()
	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	H = (np.zeros((2*numPart,2*numPart),dtype=float))#initialize Hammy with zeros, twenty by twenty in this case.

	'''More constants'''

	count = 1
	#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
	wsp_0 = mie_omegas[index]*elec/hbar + (gamma/2)*1j
	#print mie_omegas[index]
	#raw_input()
	'''initialize w_0 and eigen'''
	w_0 = 0
	eigen = np.ones(2*numPart)
	mode_count = 0
	
	
	w_mode = 0
	coupling = 0
	count = 0
	#Bfield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
	#Efield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
	for mode in [0,1]:
		vec = vectors[mode]#np.reshape(vectors[mode],(numPart,2))
		coupling = 0
		w_mode = 0
		count = 0
		while abs((w_mode*hbar/elec) - ((1 + coupling)*wsp_0*hbar/elec)) > 0.0001:
			#print abs((w_mode*hbar/elec) - ((1 + coupling)*wsp_0*hbar/elec))
			if count == 0:
				w_mode = 0
				count += 1
				#print count
			else:
				w_mode = (1 + (coupling))*wsp_0
				count += 1
				#print count
				print w_mode*hbar/elec
			wavenumber = math.sqrt(epsb)*w_mode/c
			alpha = alphasp/(1 - 1j*(2./(3.*epsb))*(wavenumber**3)*alphasp)
			coupling = 0
			near = 0
			mid = 0
			far = 0
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
					coupling += -1*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
					near += -1*r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
					mid += 1*1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
					far += -1*r_unit * exponent * (unit_dyad_term - n_dyad_term)
		w_mode = (1 + coupling)*wsp_0*hbar/elec
		if mode == 0:
			NN.append(w_mode)
			interaction[0].append(coupling)
			NF[0].append(near)
			IF[0].append(mid)
			FF[0].append(far)
		elif mode == 1:
			NS.append(w_mode)
			interaction[1].append(coupling)
			NF[1].append(near)
			IF[1].append(mid)
			FF[1].append(far)
#print len(modes[0])
#print len(modes[1])
#np.savetxt('NN_contour_water.txt',NN)
#np.savetxt('NS_contour_water.txt',NS)
r = np.linspace(1,50,50)
#dist = np.linspace(0,10,51)
plt.figure()
plt.plot(r,NN,r,NS)
plt.xlabel('radius')
plt.ylabel('energy')
plt.legend(['NN','NS'])
#plt.savefig('13_modes.pdf')
plt.show()

plt.figure()
plt.plot(r,interaction[0],label='NN',linewidth=3)
plt.plot(r,interaction[1],label='NS',linewidth=3)
plt.scatter(r,np.add(NF[0],IF[0]),label='NN near+mid',color='C0',marker='o')
plt.scatter(r,np.add(NF[1],IF[1]),label='NS near+mid',color='C1',marker='o')
plt.scatter(r,FF[0],label='NN far',color='C0',marker='s')
plt.scatter(r,FF[1],label='NS far',color='C1',marker='s')
plt.xlabel('radius')
plt.ylabel('energy')
plt.legend()
#plt.savefig('13_int.pdf')
plt.show()