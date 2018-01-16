import math
import numpy as np
import matplotlib.pyplot as plt

masses = np.loadtxt('prolate_2.5_to_1_masses.txt')
alphas = np.loadtxt('prolate_2.5_to_1_alphas.txt')
omegas = np.loadtxt('prolate_2.5_to_1_omegas.txt')
sigmas = np.loadtxt('prolate_2.5_to_1_sigmas.txt')
gammas = np.loadtxt('prolate_2.5_to_1_gammas.txt')
numPart = 12
c = 3e10
stat = 4.80326e-10
elec = 1.602e-19
hbar = 1.054e-34
epsb = 1
for r in range(1,31):
	minor = r*1e-7
	major = 2.5*r*1e-7
	mass = masses[r-1] #kg
	alpha = alphas[r-1] #m^3
	omega = omegas[r-1] #eV
	gamma = gammas[r-1] #eV
	ex = major/2
	ey = major*math.sqrt(3)/2
	Loc = [[-5*ex,0],[-4*ex,ey],[-4*ex,-ey],[-2*ex,ey],[-2*ex,-ey],[-ex,0],[ex,0],[2*ex,-ey],[2*ex,ey],[4*ex,-ey],[4*ex,ey],[5*ex,0]]
	#Loc = [[-2*major,major],[-3*major,0],[-2*major,-major],[-1*major,0],[1*major,0],[2*major,-major],[3*major,0],[2*major,major]]
	Q = [[1,0],[0,1]]
	H = (np.zeros((2*numPart,2*numPart),dtype=float))
	eigen = np.zeros(2*numPart)
	w_mode = omega
	for mode in [0,1]:
		while abs(w_mode - eigen[2*numPart-(mode+1)]) > 0.001:
			w_mode = eigen[2*numPart-(mode+1)]
			wavenumber = np.sqrt(epsb)*w_mode*elec/(c*hbar)
			for n in range(0,2*numPart):
				for m in range(0,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = (omega*elec/hbar)**2#(hbar*wsp/elec)
						#print H[n,m]
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = np.subtract(Loc[(n/2)],Loc[(m/2)]) #pick out the location of each dipole, comute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						nhat = R/Rmag #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = Rmag**-3 #this is the 1/r^3 term (static)
						r_squared = (wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (wavenumber**2)/(Rmag) #this is the 1/r term (goes with the cross products)
						space_exp = np.exp(1j*wavenumber*Rmag)
						#space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						#space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						ge = (r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))*space_exp #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -np.real(ge)*(stat**2)/(2*mass) #this has the minus sign we need.
			
			w,v = np.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting
			eigen=np.sqrt(eigenValues)*hbar/elec 
		print eigen
		
		vec = np.reshape(eigenVectors[:,2*numPart - (mode+1)],[numPart,2])
		x,y = zip(*Loc)
		u,v = zip(*vec)
		plt.figure()
		plt.title(str(r)+' nm, mode '+str(mode+1))
		plt.quiver(x,y,u,v)
		plt.show()
		#np.savetxt('hex_vec_mode_'+str(mode),vec)

