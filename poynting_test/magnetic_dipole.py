import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
from scipy.interpolate import spline
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm

a0 = 10e-8
numPart = 4
poynting_points = 32
pynts = poynting_points
mu = [0,0,1]
vec = [np.array([0,1,0])]#,np.array([1,0,0])]#,np.array([0,-1,0]),np.array([1,0,0])]
loc = [np.array([a0,0,0])]#,np.array([a0,0,0])]#,np.array([-a0,0,0]),np.array([0,-a0,0])]
screen_dist = 1000
poynting_mu = np.zeros((pynts,pynts))
'''Magnetic Dipole First'''
theta_count = 0
for theta in np.linspace(0,math.pi,pynts):
	phi_count = 0
	for phi in np.linspace(0,2*math.pi,pynts):
		nhat = [np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta)]
		point = np.multiply(screen_dist,nhat)
		rmag = np.sqrt((point[0])**2 + (point[1])**2 + point[2]**2)
		poynting_mu[theta_count,phi_count] = (np.dot(mu,mu) - (np.dot(mu,nhat))**2)
		phi_count += 1
	theta_count += 1
wavenumber = 1e5
Bfield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
Efield_total = np.zeros((poynting_points,poynting_points,3),dtype=complex)
for number in range(0,1):
	location = list(loc[number])
	#location.append(0)
	
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
			Bfield[theta_count,phi_count] = (1/rmag)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag) #
			Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)
			theta_count += 1
		phi_count += 1
	Bfield_total = Bfield_total + Bfield
	Efield_total = Efield_total + Efield
#print Bfield_total
#print Efield_total
#raw_input()
poynting_p = np.zeros((poynting_points,poynting_points),dtype=float)
for idx1 in range(poynting_points):
	for idx2 in range(poynting_points):
		poynting_p[idx1,idx2] = (rmag**2)*np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))

#print poynting_p
#raw_input()
theta = np.linspace(0,math.pi,poynting_points)
phi = np.linspace(0,2*math.pi,poynting_points)
PHI,THETA = np.meshgrid(phi,theta)
X_mu = poynting_mu * np.cos(PHI)*np.sin(THETA)
Y_mu = poynting_mu * np.sin(PHI)*np.sin(THETA)
Z_mu = poynting_mu * np.cos(THETA)
#norm = poynting_mu/poynting_mu.max()
X_p = poynting_p * np.cos(PHI)*np.sin(THETA)
Y_p = poynting_p * np.sin(PHI)*np.sin(THETA)
Z_p = poynting_p * np.cos(THETA)
norm_mu = poynting_mu/poynting_mu.max()
norm_p = poynting_p/poynting_p.max()
fig = plt.figure()
ax = fig.add_subplot(2,1,1, projection='3d')
plot = ax.plot_surface(X_mu, Y_mu, Z_mu, rstride=1, cstride=1, facecolors=cm.bwr(norm_mu),linewidth=0, antialiased=False, alpha=0.5)
#plt.axis('equal')
#plt.tight_layout()
ax = fig.add_subplot(2,1,2, projection='3d')
plot = ax.plot_surface(X_p, Y_p, Z_p, rstride=1, cstride=1, facecolors=cm.bwr(norm_p),linewidth=0, antialiased=False, alpha=0.5)
#plt.axis('equal')
#plt.tight_layout()
#plt.savefig('z_y_dipoles_test.pdf')
plt.show()