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
poynting_points = 16
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
	for r in range (20,21):
		a0 = r*10**-7; #sphere radius in cm
		index = r*10 - 10
		''' now determine geometry.'''
		print r
		#raw_input()
		# make unit vectors in centimeters.
		rij = 3*a0
		e1 = float(1)/float(2) * (rij) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (rij); #long side of 30-60-90 triangle
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
		Bfield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		Efield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
		for mode in range(0,3):
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
			screen_dist = 1000*a0
			
			for number in range(0,numPart):
				location = list(Loc[number])
				location.append(0)
				
				Bfield = np.empty((poynting_points,poynting_points,3),dtype=complex)
				Efield = np.empty((poynting_points,poynting_points,3),dtype=complex)
				phi_count = 0
				
				for phi in np.linspace(0,2*math.pi,poynting_points):
					theta_count = 0
					for theta in np.linspace(0,math.pi,poynting_points):
						nhat = [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]
						point = np.multiply(screen_dist,nhat)
						rmag = np.sqrt((point[0]-location[0])**2 + (point[1]-location[1])**2 + point[2]**2)
						nhat_dip = (location-point)/rmag
						Bfield[theta_count,phi_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)/rmag #
						Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)/math.sqrt(epsb)
						theta_count += 1
					phi_count += 1
				Bfield_total = Bfield_total + Bfield
				Efield_total = Efield_total + Efield
		poynting = np.empty((poynting_points,poynting_points),dtype=float)
		for idx1 in range(poynting_points):
			for idx2 in range(poynting_points):
				poynting[idx1,idx2] = (screen_dist**2)*np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))
		theta = np.linspace(0,math.pi,poynting_points)
		phi = np.linspace(0,2*math.pi,poynting_points)
		PHI,THETA = np.meshgrid(phi,theta)
		X = poynting * np.cos(PHI)*np.sin(THETA)
		Y = poynting * np.sin(PHI)*np.sin(THETA)
		Z = poynting * np.cos(THETA)
		norm = poynting/poynting.max()
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1, projection='3d')
		plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm),linewidth=0, antialiased=False, alpha=0.5)
		plt.show()
