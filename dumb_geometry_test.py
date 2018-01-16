import math
import numpy as np
import matplotlib.pyplot as plt

r_i = [[0,1,0],[-1,0,0],[0,-1,0],[1,0,0]]
p_i = [[-1,0,0],[0,-1,0],[1,0,0],[0,1,0]]

thing = np.zeros((181,1))

for num in range(0,4):
	theta_count = 0
	for theta in np.linspace(0,math.pi,181):
		phi_count = 0
		for phi in np.linspace(math.pi,math.pi,1):
			nhat = [np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)]
			thing[theta_count,phi_count] += np.linalg.norm(np.cross(np.cross(nhat,p_i[num]),nhat)*np.dot(nhat,r_i[num]))
			phi_count += 1
		theta_count += 1

theta = np.linspace(0,math.pi,181)
plt.figure()
plt.polar(theta,np.sin(theta))
plt.polar(theta,thing[:,0],marker='o')
plt.legend(['sin(theta)','sine and cosine term in e-field'])
#plt.savefig('sine_comparison.pdf')
plt.show()

