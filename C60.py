import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.linalg
modes = [[],[],[],[]]
mie_omegas = np.loadtxt('mie_omegas_eV.txt')
#r =1
A = np.loadtxt('buckyball_coords.txt')
for r in range(30,31):
	index = (r-1)*10
	a0 = r*1e-7
	rij = 2*a0
	ey = .5 * rij
	ex = math.sqrt(3)/2 * rij
	ez = math.sqrt(2./3.) * rij
	Q = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
	Loc = np.multiply(np.array(A[:,3:6]),rij)
	'''x,y,z = zip(*Loc)
				fig = plt.figure()
				ax = fig.gca(projection='3d')
				ax.scatter(x,y,z)
				plt.show()
				raw_input()'''
	epsb = 1
	elec = 1.60217662e-19 # regular coulombs
	numPart = len(Loc); #number of particles
	me = 9.10938291e-28; # electron mass in g
	ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
	hbar = 1.054571726e-34; # modified Planck in J*s
	c = 2.99e10; # speed of light in cm/s
	eps0 = 8.85418782e-12; # permittivity of free space
	epsinf = 3.77; # does this have units?
	Eplasma = 1.46599161e-18; # J
	gamma = 0.05*elec/(hbar*16)
	wsp_0 = mie_omegas[index] * elec/hbar
	alpha_stat = (3*a0**3)/(epsinf+2*epsb)
	w_0 = 0
	eigen = np.ones(3*numPart)
	H = np.zeros((3*numPart,3*numPart))
	count = 1
	for mode in range(0,1):
		while np.sqrt(np.square(w_0*hbar/elec - eigen[(3*numPart)-(mode+1)])) > 0.0000001:
			if count == 1:
				wsp = wsp_0
				count = count + 1
				w_0 = 0
			else:
				count = count + 1
				w_0 = eigen[3*numPart-(mode+1)]*elec/hbar
			wavenumber = (w_0*math.sqrt(epsb))/(c)
			alpha = alpha_stat/(1 - 1j*(2./3.)*(wavenumber**3)*alpha_stat)
			msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
			gamma_eV = gamma_ret*hbar/elec
			#wsp = wsp_0
			wsp = math.sqrt((np.real(wsp_0))**2 - (np.real(gamma_ret)/2)**2) # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,3*numPart):
				for m in range (0,3*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1#(hbar/elec)*wsp
					elif np.sqrt(np.sum(np.square(Loc[(n/3)]-Loc[(m/3)]))) < a0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/3)]-Loc[(m/3)] #pick out the location of each dipole, compute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2+R[2]**2) #compute magnitude of distance
						nhat = (Loc[(n/3)]-Loc[(m/3)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%3],Q[m%3]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%3],nhat)*np.dot(nhat,Q[m%3]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						exponent = np.exp(1j*w_0*Rmag/c)
						ge = ((r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))) * exponent #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -np.real(ge)/epsb#*(hbar/elec)*wsp #this has the minus sign we need.
						#H[m,n] = -np.real(ge)/epsb
						#H[m,n] = np.conj(-ge)
			#print H
			#raw_input()
			diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
			Ht = np.matrix.transpose(H) # this is the transpose of H
			Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			#print Hfull
			#raw_input()
			w,v = scipy.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			#print idx
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting

			eigen=((hbar/elec)*wsp_0)*(np.sqrt(eigenValues))# the eigenvalues have units of energy^2, so we take the square root
		modes[mode].append(eigen[3*numPart-(mode+1)])
		vec = np.reshape(eigenVectors[:,3*numPart-(mode+1)],(numPart,3))
		
print eigen[-20:-1]
#raw_input()
vec = np.reshape(eigenVectors[:,-9],(numPart,3))
x,y,z = zip(*Loc)
u,v,w = zip(*vec)
fig = plt.figure()
#plt.scatter(x,y)
#plt.quiver(x,y,u,v)
ax = fig.gca(projection='3d')
ax.scatter(x,y,z)
ax.quiver(x,y,z,u,v,w,length=5*rij)
plt.axis('scaled')
plt.show()
#sep = np.linspace(1,30,30)
#plt.figure()
#plt.plot(sep,modes[0],sep,modes[1],sep,modes[2],sep,modes[3],linewidth=3)
#plt.savefig('separating_'+str(r)+'_nm_rij=2.5a0.pdf')
#plt.show()



'''while np.sqrt(np.square(w_0*hbar/elec - eigen[(2*numPart)-(mode+1)])) > 0.0000001:
			if count == 1:
				wsp = wsp_0
				count = count + 1
				w_0 = 0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
			wavenumber = (w_0*math.sqrt(epsb))/(c)
			alpha = alpha_stat/(1 - 1j*(2./3.)*(wavenumber**3)*alpha_stat)
			msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping time
			gamma_ret = gamma+tau*(w_0**2) # I pulled this from the beats paper
			gamma_eV = gamma_ret*hbar/elec
			#wsp = wsp_0
			wsp = math.sqrt((np.real(wsp_0))**2 - (np.real(gamma_ret)/2)**2) # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,3*numPart):
				for m in range (0,3*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1#(hbar/elec)*wsp
					elif np.sqrt(np.sum(np.square(Loc[(n/3)]-Loc[(m/3)]))) < a0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/3)]-Loc[(m/3)] #pick out the location of each dipole, compute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2+R[2]**2) #compute magnitude of distance
						nhat = (Loc[(n/3)]-Loc[(m/3)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%3],Q[m%3]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%3],nhat)*np.dot(nhat,Q[m%3]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*w_0)/(c*(Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*w_0**2)/(Rmag*(c**2)) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						exponent = np.exp(1j*w_0*Rmag/c)
						ge = ((r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))) * exponent #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -np.real(ge)/epsb#*(hbar/elec)*wsp #this has the minus sign we need.
						#H[m,n] = -np.real(ge)/epsb
						#H[m,n] = np.conj(-ge)
			#print H
			#raw_input()
			diag = np.diag(np.diag(H)) # this produces a matrix of only the diagonal terms of H
			Ht = np.matrix.transpose(H) # this is the transpose of H
			Hedit = diag - Ht # this produces a matrix with zeros on the diagonal and the upper triangle, and the lower triangle has all the leftover values of H with the opposite sign
			Hfull = H - Hedit # this combines H with the lower triangle (all negative) to produce a symmetric, full matrix
			#print Hfull
			#raw_input()
			w,v = scipy.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			#print idx
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting

			eigen=((hbar/elec)*wsp_0)*(np.sqrt(eigenValues))# the eigenvalues have units of energy^2, so we take the square root
		modes[mode].append(eigen[3*numPart-(mode+1)])
		vec = np.reshape(eigenVectors[:,3*numPart-(mode+1)],(numPart,3))
		np.savetxt('bitetra_eigen_' + str(mode),vec)'''