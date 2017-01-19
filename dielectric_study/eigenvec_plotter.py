import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt

for r in range(1,21):
	numPart = 10; #number of particles
	a0 = r*10**-7; #sphere radius in cm

	e1 = float(1)/float(2) * (2.2*a0) ; #short side of 30-60-90 triangle
	e2 = float(math.sqrt(3))/float(2) * (2.2*a0); #long side of 30-60-90 triangle
	Loc = [np.array([0, e1]),np.array([-e2, 2*e1]),np.array([-2*e2, e1]),np.array([-2*e2, -e1]),np.array([-e2, -2*e1]),
		   np.array([0, -e1]),np.array([e2, -2*e1]),np.array([2*e2, -e1]),np.array([2*e2 , e1]),np.array([e2 , 2*e1])] #location vectors for center of each sphere - they go counterclockwise, with one being the top center particle and six being the bottom center particle.


	eigenvec = np.loadtxt('/'.join([str(1.42),'_'.join([str(r),'nm_lowest.txt'])]))
	vec = np.reshape(eigenvec, (10,2))
	#print Loc[:]
	x,y = zip(*Loc)
	#y = Loc[:,1]
	u,v = zip(*vec)
	#v = vec[:,1]
	plt.title(''.join(['radius =',str(r)]))
	plt.quiver(x,y,u,v)
	plt.show()