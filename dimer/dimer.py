import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
coll_near=[]
coll_int=[]
coll_far =[]
coll_tot = []

anpa_near = []
anpa_int = []
anpa_far = []
anpa_tot = []

para_near = []
para_int = []
para_far = []
para_tot = []

anco_near = []
anco_int =[]
anco_far = []
anco_tot = []
collinear = [np.array([1,0]),np.array([1,0])]
antiparallel = [np.array([0,1]),np.array([0,-1])]
parallel = [np.array([0,1]),np.array([0,1])]
anticollinear = [np.array([1,0]),np.array([-1,0])]
r = 10
elec = 1.60217662e-19 # regular coulombs
numPart = 2; #number of particles
a0 = 30e-7; #sphere radius in cm

''' now determine geometry.'''

# make unit vectors in centimeters.
for number in range(0,481):
	#dist = 
	sep = ((number+20)*1e-7)+2*a0
	Loc = [np.array([sep/2, 0]),np.array([-sep/2, 0])]

	#Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction

	'''This part builds the Hamiltonian. We start by initializing and choosing an initial omega value.'''
	#H = np.zeros((2*numPart,2*numPart)) #initialize Hammy with zeros, twenty by twenty in this case.

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
	wsp_0 = mie_omegas[(r-1)*10]*elec/hbar
	'''initialize w_0 and eigen'''
	w_0 = wsp_0
	eigen = np.ones(2*numPart)
	wavenumber = (w_0)/(c*math.sqrt(epsb))
	alphasp = (a0**3)*(3/(epsinf+2*epsb))
	alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
	R = Loc[(1)]-Loc[(0)] #pick out the location of each dipole, comute the distance between them
	Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
	nhat = (Loc[(1)]-Loc[(0)])/float(Rmag)
	r_cubed = alpha/Rmag**3 #this is the 1/r^3 term (static)
	r_squared = (w_0*alpha)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
	r_unit = (alpha*w_0**2)/(Rmag*(c**2))
	space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
	space_sin = np.sin(w_0*Rmag/c)

	coll_near.append(-(3*np.dot(collinear[0],nhat)*np.dot(collinear[1],nhat) - np.dot(collinear[0],collinear[1]) ) * r_cubed * space_cos)
	anpa_near.append(-(3*np.dot(antiparallel[0],nhat)*np.dot(antiparallel[1],nhat) - np.dot(antiparallel[0],antiparallel[1]) ) * r_cubed * space_cos)
	anco_near.append(-(3*np.dot(anticollinear[0],nhat)*np.dot(anticollinear[1],nhat) - np.dot(anticollinear[0],anticollinear[1]) ) * r_cubed * space_cos)
	para_near.append(-(3*np.dot(parallel[0],nhat)*np.dot(parallel[1],nhat) - np.dot(parallel[0],parallel[1]) ) * r_cubed * space_cos)

	coll_int.append(-(3*np.dot(collinear[0],nhat)*np.dot(collinear[1],nhat) - np.dot(collinear[0],collinear[1]) ) * r_squared * space_sin)
	anpa_int.append(-(3*np.dot(antiparallel[0],nhat)*np.dot(antiparallel[1],nhat) - np.dot(antiparallel[0],antiparallel[1]) ) * r_squared * space_sin)
	anco_int.append(-(3*np.dot(anticollinear[0],nhat)*np.dot(anticollinear[1],nhat) - np.dot(anticollinear[0],anticollinear[1]) ) * r_squared * space_sin)
	para_int.append(-(3*np.dot(parallel[0],nhat)*np.dot(parallel[1],nhat) - np.dot(parallel[0],parallel[1]) ) * r_squared * space_sin)

	coll_far.append((np.dot(collinear[0],nhat)*np.dot(collinear[1],nhat) - np.dot(collinear[0],collinear[1]) ) * r_unit * space_cos)
	anpa_far.append((np.dot(antiparallel[0],nhat)*np.dot(antiparallel[1],nhat) - np.dot(antiparallel[0],antiparallel[1]) ) * r_unit * space_cos)
	anco_far.append((np.dot(anticollinear[0],nhat)*np.dot(anticollinear[1],nhat) - np.dot(anticollinear[0],anticollinear[1]) ) * r_unit * space_cos)
	para_far.append((np.dot(parallel[0],nhat)*np.dot(parallel[1],nhat) - np.dot(parallel[0],parallel[1]) ) * r_unit * space_cos)

	coll_tot.append(coll_near[number] + coll_int[number] + coll_far[number])
	anpa_tot.append(anpa_near[number] + anpa_int[number] + anpa_far[number])
	anco_tot.append(anco_near[number] + anco_int[number] + anco_far[number])
	para_tot.append(para_near[number] + para_int[number] + para_far[number])
#print coll_tot
sep = np.linspace(20,500,481)
plt.figure(1)
plt.plot(sep,coll_near,sep,coll_int,sep,coll_far,sep,coll_tot,lw=3)
plt.xlabel('Separation Distance (nm)')
plt.ylabel('Interaction Strength')
#plt.xlim([2.5,250])
plt.legend(['near','int','far','total'],loc=4)
#plt.savefig('collinear.pdf')
plt.show()

plt.figure(2)
plt.plot(sep,anpa_near,sep,anpa_int,sep,anpa_far,sep,anpa_tot,lw=3)
plt.xlabel('Separation Distance (nm)')
plt.ylabel('Interaction Strength')
#plt.xlim([2.5,250])
plt.legend(['near','int','far','total'],loc=4)
#plt.savefig('antiparallel.pdf')
plt.show()

plt.figure(3)
plt.plot(sep,anco_near,sep,anco_int,sep,anco_far,sep,anco_tot,lw=3)
plt.xlabel('Separation Distance (nm)')
plt.ylabel('Interaction Strength')
#plt.xlim([2.5,250])
plt.legend(['near','int','far','total'])
#plt.savefig('anticollinear.pdf')
plt.show()

plt.figure(4)
plt.plot(sep,para_near,sep,para_int,sep,para_far,sep,para_tot,lw=3)
plt.xlabel('Separation Distance (nm)')
plt.ylabel('Interaction Strength')
#plt.xlim([2.5,250])
plt.legend(['near','int','far','total'])
#plt.savefig('parallel.pdf')
plt.show()

#plt.show()

