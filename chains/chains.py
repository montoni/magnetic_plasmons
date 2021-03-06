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
	for r in range (1,2):
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

		count = 1
		#wsp_0 = math.sqrt((wplasma/math.sqrt(epsinf+2*epsb))**2 - (gamma/2)**2);
		wsp_0 = (mie_omegas[index])*elec/hbar
		'''initialize w_0 and eigen'''
		w_0 = wsp_0
		eigen = np.zeros(2*numPart)
		for mode in range(2,3):
			mag_dipole = []
			while abs(w_0*hbar/elec - eigen[(2*numPart)-(mode+1)]) > 0.00001:
				if count == 1:
					wsp = wsp_0
					count = count + 1
					w_0 = 0
				else:
					count = count + 1
					w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
				#print w_0
				alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
				wavenumber = (w_0)/(c*math.sqrt(epsb))
				alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
				msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
				tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
				gamma_ret = gamma+tau*(wsp**2) # I pulled this from the beats paper
				gamma_eV = gamma_ret*hbar/elec
				#wsp = wsp_0
				wsp = wsp_0 # sp frequency (rad/s) corrected for radiation damping
				for n in range (0,2*numPart):
					for m in range (0,2*numPart):
						if m == n: #if m and n are the same, the hammy gets the plasmon energy
							H[n,m] = 1
						elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
							H[n,m] = 0
						elif n == m+1 and m%2 == 0:
							H[m,n] = 0
						else: # all other dipoles are fair game
							R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, comute the distance between them
							Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
							nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
							p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
							p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
							r_cubed = alpha/Rmag**3 #this is the 1/r^3 term (static)
							r_squared = (-1j*alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
							r_unit = (alpha*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
							space_exp = np.exp(1j*w_0*Rmag/c)
							#space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
							#space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
							ge = (r_unit * (p_dot_p - p_nn_p) + (r_cubed + r_squared) * (3*p_nn_p - p_dot_p))*space_exp #this is p dot E
							#print ge
							gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
							H[n,m] = - np.real(ge) #this has the minus sign we need.
							H[m,n] = - np.real(ge)
				'''diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
				Ht = np.matrix.transpose(H) # this is the transpose of H
				Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
				Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix'''
				w,v = np.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
				idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
				eigenValues = w[idx] # sorting
				eigenVectors = v[:,idx] # sorting
				eigen=np.sqrt(eigenValues)*wsp*hbar/elec # the eigenvalues have units of energy^2, so we take the square root
				#print eigen
			    #w_old = w_0
			    #w_0 = eigen[2*numPart-1]
			vec = np.reshape(eigenVectors[:,2*numPart - (mode+1)],[numPart,2])
			grid_x = 4.5*e2
			grid_y = 2.5*e1
			size = 251
			Bfield_total = np.zeros((size,size),dtype=float)
			for number in range(0,numPart):
				point = Loc[number]
				xcount = 0
				Bfield = np.empty((size,size),dtype=float)
				for xx in np.linspace(-grid_y,grid_y,size):
					ycount = 0
					for yy in np.linspace(-grid_x,2.5*e2,size):
						rmag = np.sqrt((point[0]-yy)**2 + (point[1]-xx)**2)
						#print rmag
						#raw_input()
						#vec[number][0],vec[number][1] = vec[number][1],vec[number][0]
						if rmag < .5*a0:
							Bfield[xcount,ycount] = 0
						else:
							nhat = -([yy,xx]-point)/rmag
							Bfield[xcount,ycount] = (wavenumber**2)*np.cross(nhat,vec[number])*(rmag**-1)*(1-((1j*wavenumber*rmag)**-1))#*np.exp(1j*wavenumber*rmag)
						ycount += 1
					xcount += 1
				Bfield_total = Bfield_total + Bfield
			# ,levels=np.linspace(np.amin(Bfield_total),np.amax(Bfield_total),30)
			xx = np.linspace(-grid_x,2.5*e2,size)
			yy = np.linspace(-grid_y,grid_y,size)
			#Loc = places
			u,v = zip(*vec)
			xloc, yloc = zip(*Loc)
			plt.figure()
			plt.contourf(xx,yy,Bfield_total,cmap='bwr',levels=np.linspace(np.amin(Bfield_total),np.amax(Bfield_total),30))
			#plt.colorbar()
			plt.tight_layout()
			plt.scatter(xloc,yloc)
			plt.quiver(xloc,yloc,u,v)
			plt.savefig('mode_'+str(mode)+ '_field.pdf')
			plt.show()
			'''screen_dist = 20*a0
			Bfield_total = np.zeros((361,3),dtype=complex)
			Efield_total = np.zeros((361,3),dtype=complex)
			for number in range(0,numPart):
				location = list(Loc[number])
				location.append(0)
				
				Bfield = np.empty((361,3),dtype=complex)
				Efield = np.empty((361,3),dtype=complex)
				theta_count = 0
				for theta in np.linspace(0,2*math.pi,361):
					nhat = [0, np.cos(theta), np.sin(theta)]
					point = np.multiply(screen_dist,nhat)
					
					rmag = np.sqrt((point[0]-location[0])**2 + (point[1]-location[1])**2 + point[2]**2)
					nhat_dip = (location-point)/rmag
					Bfield[theta_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)
					Efield[theta_count] = np.cross(Bfield[theta_count],nhat_dip)
					theta_count += 1
				Bfield_total = Bfield_total + Bfield
				Efield_total = Efield_total + Efield
			poynting = np.empty((361,1),dtype=complex)
			for idx in range(361):
				poynting[idx] = np.linalg.norm(np.real(np.cross(Efield_total[idx],np.conj(Bfield_total[idx]))))
			theta = np.linspace(0,2*math.pi,361)
			plt.figure()
			plt.polar(theta,poynting)
			plt.show()'''