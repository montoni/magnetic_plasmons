import math
import numpy as np
import matplotlib.pyplot as plt

masses = np.loadtxt('prolate_2.5_to_1_masses.txt')
alphas = np.loadtxt('prolate_2.5_to_1_alphas.txt')
omegas = np.loadtxt('prolate_2.5_to_1_omegas.txt')
sigmas = np.loadtxt('prolate_2.5_to_1_sigmas.txt')
gammas = np.loadtxt('prolate_2.5_to_1_gammas.txt')
numPart = 6
c = 3e10
stat = 4.80326e-10
elec = 1.602e-19
hbar = 1.054e-34
epsb = 1
NN = []
NS = []
interaction = [[],[]]
NF = [[],[]]
IF = [[],[]]
FF = [[],[]]
NN_vec = np.loadtxt('vectors_mode_1')
NS_vec = np.loadtxt('vectors_mode_2')
vector = [NN_vec,NS_vec]
for r in range(1,31):
	minor = r*1e-7
	major = 2.5*r*1e-7
	mass = masses[r-1] #kg
	alpha = alphas[r-1] #m^3
	omega = omegas[r-1] #eV
	wsp = omega*elec/hbar
	gamma = gammas[r-1] #eV
	ex = (math.sqrt(3) + 2)*major
	ey = major
	Loc = [[-ex,ey],[-ex,-ey],[-2*major,0],[2*major,0],[ex,-ey],[ex,ey]]
	#Loc = [[-2*major,major],[-3*major,0],[-2*major,-major],[-1*major,0],[1*major,0],[2*major,-major],[3*major,0],[2*major,major]]
	Q = [[1,0],[0,1]]
	w_mode = 0
	coupling = 0
	for mode in [1,2]:
		vec = vector[mode-1]
		w_mode = 0
		coupling = 0
		while abs(w_mode - (omega + np.real(coupling))) > 0.001:
			w_mode = omega + np.real(coupling)
			wavenumber = np.sqrt(epsb)*w_mode*elec/(c*hbar)
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
					R = np.subtract(Loc[(x)],Loc[(y)]) #pick out the location of each dipole, compute the distance between them
					Rmag = math.sqrt(R[0]**2+R[1]**2) #compute magnitude of distance
					unit_vector = R/Rmag #compute unit vector between dipoles
					unit_dyad_term = epsb*np.dot(vec[x],vec[y])
					n_dyad_term = epsb*np.dot(vec[x],unit_vector)*np.dot(unit_vector,vec[y])
					r_cubed = 1/(Rmag**3) #this is the 1/r^3 term (static)
					r_squared = (wavenumber)/((Rmag**2)) #this is the 1/r^2 term (imaginary part)
					r_unit = (wavenumber**2)/(Rmag)
					exponent = np.exp(1j*wavenumber*Rmag)
					coupling += -(hbar*(stat)**2)/(2*mass*wsp*elec)*(((r_unit * (unit_dyad_term - n_dyad_term) + (r_cubed - 1j*r_squared) * (3*n_dyad_term - unit_dyad_term))) * exponent)
					nearfield += -(hbar*(stat)**2)/(mass*wsp*elec)*r_cubed*exponent*(3*n_dyad_term - unit_dyad_term)
					midfield += (hbar*(stat)**2)/(mass*wsp*elec)*1j*r_squared * (3*n_dyad_term - unit_dyad_term) * exponent
					farfield += -(hbar*(stat)**2)/(mass*wsp*elec)*r_unit * exponent * (unit_dyad_term - n_dyad_term)
			
		if mode == 1:
			NN.append(w_mode)
		else:
			NS.append(w_mode)
print NN
print NS
r = np.linspace(1,30,30)
plt.figure()
plt.plot(r,NN,r,NS)
plt.legend(['NN','NS'])
plt.show()
		
		

