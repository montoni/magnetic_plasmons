import math
import numpy as np

omega_sp = [2.85,2.84,2.84,2.83,2.83,2.82,2.81,2.8,2.78,2.77,2.75,2.74,2.72,2.7,2.68,2.66,2.64,2.61,2.59,2.56,2.54,2.51,2.48,2.46,2.43,2.40,2.37,2.34,2.32,2.29]
sigma_abs = [44.76,352.2,1202,2670,5051,8126,11610,15080,18340,21350,23220,24540,25250,25280,24830,24050,23030,22080,21000,19910,18880,17910,16960,16150,15380,14680,14040,13460,12960,12500]
gamma_sp = [0.05,0.05,0.05,0.05,0.05,0.05,0.06,0.06,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.16,0.18,0.19,0.2,0.22,0.24,0.26,0.27,0.3,0.33,0.34,0.37,0.4]
msp = []
alpha_sp = []
stat = 4.80326e-10
coul = 1.602e-19
hbar = 1.054e-34
c = 3e10
for r in range(1,31):
	mass = (2*math.pi*stat**2)/(3*c*sigma_abs[r-1]*(1e-14)*gamma_sp[r-1]*coul/hbar)
	alpha = (stat**2)/(mass*(omega_sp[r-1]**2)*(coul/hbar)**2)
	msp.append(mass)
	alpha_sp.append(alpha)

np.savetxt('prolate_2.5_to_1_masses.txt',msp)
np.savetxt('prolate_2.5_to_1_alphas.txt',alpha_sp)
np.savetxt('prolate_2.5_to_1_omegas.txt',omega_sp)
np.savetxt('prolate_2.5_to_1_sigmas.txt',sigma_abs)
np.savetxt('prolate_2.5_to_1_gammas.txt',gamma_sp)