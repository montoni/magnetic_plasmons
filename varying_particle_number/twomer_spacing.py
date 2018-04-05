import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib.pyplot import cm

mie_omegas = np.loadtxt('../mie_omegas_eV.txt')
#mie_omegas = np.loadtxt('../mie_frequencies_silver_water.txt')
''' So all the units are cgs, but the hamiltonian gets loaded up with energies in eVs, so the first constant below is the charge of an electron in coulombs and the rest of the units are cgs. Every constant should be labeled.'''
crossing = []
elec = 1.60217662e-19 # regular coulombs
maxlength = []
numPart = 10; #number of particles
NS_vec = np.loadtxt('ten_mode_0.txt')
NN_vec = np.loadtxt('ten_mode_1.txt')
vectors = [[NS_vec],[NN_vec]]
me = 9.10938291e-28; # electron mass in g
ch = 4.80326e-10; # electron charge in g^1/2 cm^3/2 s^-1
hbar = 1.054571726e-34; # modified Planck in J*s
c = 2.99e10; # speed of light in cm/s
eps0 = 8.85418782e-12; # permittivity of free space
Q = [np.array([1,0]),np.array([0,1])] #dipole moments: 1st in x-direction, 2nd in y-direction
epsinf = 3.77; 
'''Properties for silver.'''
Eplasma = 1.46599161e-18; # J
gamma = 0.05*elec/(hbar*16)
wplasma = Eplasma/hbar; # plasma frequency (rad/s)
epsb = 1
NN = np.zeros((50,51))
NS = np.zeros((50,51))
modes = [[],[]]
interaction = [[],[]]
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
dist = 1
r_count = 0
for r in np.linspace(1,50,51):
	r_count += 1
	dist_count = 0
	print r
#r = 15
	for dist in np.linspace(1,1,1):
		dist_count += 1

		a0 = r*10**-7; #sphere radius in cm
		alphasp = (a0**3)*(3/(epsinf+2*epsb)); # polarizability (cm^3)
		index = r
		
		''' now determine geometry.'''
		print dist
		# make unit vectors in centimeters.
		rij = (dist+2)*a0
		part_per_ring = numPart/2 + 1
		theta = 2*math.pi/part_per_ring
		phi = theta/2.
		#print theta
		#print phi
		#print np.cos(phi)
		#print np.sin(phi)
		cent_dist = rij/(2*np.tan(phi))
		part_to_cent = math.sqrt((rij/2)**2 + (cent_dist)**2)
		centers = [np.array([-cent_dist,0]),np.array([cent_dist,0])]
		#print centers
		places = []
		for num in range(part_per_ring-1):
			if part_per_ring%2 == 0:
				#print 'hello'
				places.append(centers[0] + np.array([(part_to_cent*np.cos(phi+theta*(num))),(part_to_cent*np.sin(phi+theta*(num)))]))
				places.append(centers[1] + np.array([(part_to_cent*np.cos(phi+theta*(num+np.ceil(part_per_ring/2.)))),(part_to_cent*np.sin(phi+theta*(num+np.ceil(part_per_ring/2.))))]))
			elif part_per_ring%2 == 1:
				#print 'good bye'
				places.append(centers[0] + np.array([(part_to_cent*np.cos(phi+theta*(num))),(part_to_cent*np.sin(phi+theta*(num)))]))
				places.append(centers[1] + np.array([(part_to_cent*np.cos(theta*(num+np.ceil(part_per_ring/2.)))),(part_to_cent*np.sin(theta*(num+np.ceil(part_per_ring/2.))))]))
		'''Loc = np.zeros((numPart,2))
		count = 0 
		for x in places:
			print x
			if np.isclose(x.all(),Loc[:].any()):
				pass
			else:
				Loc[count] = x'''


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
		wsp_0 = mie_omegas[r_count-1] * elec/hbar
		#wsp_0 = (mie_omegas[r-1]-0.05)*elec/hbar -r*0.0034
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
					#print w_mode*hbar/elec
				wavenumber = math.sqrt(epsb)*w_mode/c
				alpha = ((alphasp**-1) - 1j*(2./(3.*epsb))*wavenumber**3)**-1
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
						coupling += -(hbar/elec)*wsp_0*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
						near += -(hbar/elec) * wsp_0 * r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
						mid += -(hbar/elec) * wsp_0 * -1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
						far += -(hbar/elec) * wsp_0 * r_unit * exponent * (unit_dyad_term - n_dyad_term)
			w_mode = wsp_0*hbar/elec + np.real(coupling)
			if mode == 0:
				NS[r_count-1,dist_count-1] = w_mode
				'''interaction[1].append(coupling)
				NF[1].append(near)
				IF[1].append(mid)
				FF[1].append(far)'''
			elif mode == 1:
				NN[r_count-1,dist_count-1] = w_mode
				'''interaction[0].append(coupling)
				NF[0].append(near)
				IF[0].append(mid)
				FF[0].append(far)'''
#print len(modes[0])
#print len(modes[1])
np.savetxt('NN_contour_water.txt',NN)
np.savetxt('NS_contour_water.txt',NS)
r = np.linspace(1,30,291)
dist = np.linspace(0,10,51)
plt.figure()
plt.plot(dist,NN[0],dist,NS[0])
plt.legend(['NN','NS'])
plt.show()
#np.savetxt('15nm_NN_spacing_water.txt',np.real(modes[0]))
#np.savetxt('15nm_NS_spacing_water.txt',np.real(modes[1]))
#del(modes[0][50])
#plt.figure()
#plt.contourf(NN-NS)
#plt.colorbar()
#plt.show()
'''plt.plot(dist,NN)
plt.plot(dist,NS)
plt.plot(dist,np.multiply(maxlength,1e7))
plt.ylabel('Wavelength (nm)')
plt.xlabel('s/r_0')		
plt.legend(['NN','NS'])
plt.savefig('water_twomer_scale_eig.pdf')
plt.show()'''


'''plt.figure()
plt.plot(dist,interaction[0],label='NN',linewidth=3)
plt.plot(dist,interaction[1],label='NS',linewidth=3)
plt.scatter(dist,np.add(NF[0],IF[0]),label='NN near+mid',color='C0',marker='o')
plt.scatter(dist,np.add(NF[1],IF[1]),label='NS near+mid',color='C1',marker='o')
plt.scatter(dist,FF[0],label='NN far',color='C0',marker='s')
plt.scatter(dist,FF[1],label='NS far',color='C1',marker='s')
plt.ylabel('Energy (eV)')
plt.xlabel('Particle Radius r_0 (nm)')	
plt.legend()
#plt.savefig('water_twomer_scale_int.pdf')
plt.show()'''
#plt.savefig('twomer_interactions.pdf')

'''screen_dist = 1000*a0
			
			for number in range(0,numPart):
				location = list(Loc[number])
				location.append(0)
				
				Bfield = np.empty((16,16,3),dtype=complex)
				Efield = np.empty((16,16,3),dtype=complex)
				phi_count = 0
				
				for phi in np.linspace(0,2*math.pi,16):
					theta_count = 0
					for theta in np.linspace(0,math.pi,16):
						nhat = [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]
						point = np.multiply(screen_dist,nhat)
						rmag = np.sqrt((point[0]-location[0])**2 + (point[1]-location[1])**2 + point[2]**2)
						nhat_dip = (location-point)/rmag
						Bfield[theta_count,phi_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)
						Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)
						theta_count += 1
					phi_count += 1
				Bfield_total = Bfield_total + Bfield
				Efield_total = Efield_total + Efield
		poynting = np.empty((16,16),dtype=float)
		for idx1 in range(16):
			for idx2 in range(16):
				poynting[idx1,idx2] = np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))
		theta = np.linspace(0,math.pi,16)
		phi = np.linspace(0,2*math.pi,16)
		PHI,THETA = np.meshgrid(phi,theta)
		X = poynting * np.cos(PHI)*np.sin(THETA)
		Y = poynting * np.sin(PHI)*np.sin(THETA)
		Z = poynting * np.cos(THETA)
		norm = poynting/poynting.max()
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1, projection='3d')
		plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm),linewidth=0, antialiased=False, alpha=0.5)
		plt.show()'''


'''for mode in range(0,8):
		mag_dipole = []

		while np.sqrt(np.square(w_0*hbar/elec - eigen[(2*numPart)-(mode+1)])) > 0.00001:
			if count == 1:
				wsp = wsp_0
				count = count + 1
				w_0 = 0
			else:
				count = count + 1
				w_0 = eigen[2*numPart-(mode+1)]*elec/hbar
			wavenumber = (w_0*math.sqrt(epsb))/(c)
			alpha = alphasp/(1 - 1j*(2./3.)*(wavenumber**3)*alphasp)
			msp = (ch**2)/(alpha*((wsp)**2)); # sp mass (grams)
			tau = (2*ch**2)/(3*msp*c**3) # radiation damping period
			gamma_ret = gamma+tau*(w_0**2) 
			gamma_eV = gamma_ret*hbar/elec
			#wsp = wsp_0
			#wsp = math.sqrt((wsp_0)**2 - (gamma_ret/2)**2) # sp frequency (rad/s) corrected for radiation damping
			for n in range (0,2*numPart):
				for m in range (0,2*numPart):
					if m == n: #if m and n are the same, the hammy gets the plasmon energy
						H[n,m] = 1#(hbar/elec)*wsp
					elif m == n+1 and n%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					elif n == m+1 and m%2 == 0: #if m and n are on the same particle, they don't couple
						H[n,m] = 0
					else: # all other dipoles are fair game
						R = Loc[(n/2)]-Loc[(m/2)] #pick out the location of each dipole, compute the distance between them
						Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
						distances.append(Rmag)
						nhat = (Loc[(n/2)]-Loc[(m/2)])/float(Rmag) #compute unit vector between dipoles
						p_dot_p = np.dot(Q[n%2],Q[m%2]) # this is one dipole dotted into the other
						p_nn_p = np.dot(Q[n%2],nhat)*np.dot(nhat,Q[m%2]) # this is one dipole dotted into the unit vector, times the other dipole dotted into the unit vector
						r_cubed = alpha/(Rmag**3) #this is the 1/r^3 term (static)
						r_squared = (alpha*wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
						r_unit = (alpha*wavenumber**2)/(Rmag) #this is the 1/r term (goes with the cross products)
						#space_exp = np.exp(1j*w_0*Rmag/c)
						space_cos = np.cos(w_0*Rmag/c) #this is the real part of the e^ikr
						space_sin = np.sin(w_0*Rmag/c) #this is the imaginary part of the e^ikr
						exponent = np.exp(1j*wavenumber*Rmag)
						ge = ((r_unit * (p_dot_p - p_nn_p) + (r_cubed - 1j*r_squared) * (3*p_nn_p - p_dot_p))) * exponent #this is p dot E
						gm = 0 #set magnetic coupling to zero. we can include this later if necessary.
						H[n,m] = -np.real(ge)/epsb#*(hbar/elec)*wsp #this has the minus sign we need.
						
			w,v = scipy.linalg.eig(H) #this solves the eigenvalue problem, producing eigenvalues w and eigenvectors v.
			idx = w.argsort()[::-1] # this is the idx that sorts the eigenvalues from largest to smallest
			#print idx
			eigenValues = w[idx] # sorting
			eigenVectors = v[:,idx] # sorting

			eigen=((hbar/elec)*wsp)*(np.sqrt(eigenValues))# the eigenvalues have units of energy^2, so we take the square root
			vec = np.reshape(eigenVectors[:,2*numPart - (mode+1)],[numPart,2])'''
