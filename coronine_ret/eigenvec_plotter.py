import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

numPart = 24; #number of particles

for r in range(18,41):
	for num in range(1,8):
		a0 = 0.5*r*10**-7; #sphere radius in cm

		e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
		e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
		Loc = [np.array([0,2*e1]),np.array([-e2,e1]),np.array([-e2,-e1]),np.array([0,-2*e1]),
			   np.array([e2,-e1]),np.array([e2,e1]),np.array([2*e2,2*e1]),np.array([3*e2,e1]),
			   np.array([3*e2,-e1]),np.array([2*e2,-2*e1]),np.array([2*e2,-4*e1]),np.array([e2,-5*e1]),
			   np.array([0,-4*e1]),np.array([-e2,-5*e1]),np.array([-2*e2,-4*e1]),np.array([-2*e2,-2*e1]),
			   np.array([-3*e2,-e1]),np.array([-3*e2,e1]),np.array([-2*e2,2*e1]),np.array([-2*e2,4*e1]),
			   np.array([-e2,5*e1]),np.array([0, 4*e1]),np.array([e2,5*e1]),np.array([2*e2,4*e1])]

		x,y = zip(*Loc)
		eigenvec = np.loadtxt('_'.join([str(r),'nm','.'.join([str(num),'txt'])]))
		vec = np.reshape(eigenvec, (24,2))
		#print Loc[:]
		
		#y = Loc[:,1]
		u,v = zip(*vec)
		#v = vec[:,1]
		plt.title(''.join(['radius =',str(r*0.5)]))
		plt.quiver(x,y,u,v)
		plt.show()