#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
'''NN_eig = np.loadtxt('NN_eig')
NS_eig = np.loadtxt('NS_eig')
r = [1,2.5,5,6,7,8,9,10]
plt.plot(r,NN_eig,r,NS_eig,linewidth=3)
plt.show()'''


'''cross = np.loadtxt('crossing')
eps = np.arange(1,16)
plt.plot(eps,np.square(cross),linewidth=3)
plt.show()'''
NN_surface = np.zeros([10,291])
NS_surface = np.zeros([10,291])
X = []
Y = []
for epsb in range (1,11):
	A = np.loadtxt('.'.join(['_'.join(['NN_eps',str(epsb)]),'txt']))
	B = np.loadtxt('.'.join(['_'.join(['NS_eps',str(epsb)]),'txt']))
	NN_surface[epsb-1] = A
	NS_surface[epsb-1] = B

x = np.linspace(1,30,291)
epsb = np.linspace(1,10,10)
limit = np.linspace(1,10,10)

diff_5 = (NN_surface[:,51]-NS_surface[:,51])
diff_10 = (NN_surface[:,101]-NS_surface[:,101])
diff_15 = (NN_surface[:,151]-NS_surface[:,151])
diff_20 = (NN_surface[:,201]-NS_surface[:,201])

qnew = np.linspace(epsb.min(),epsb.max(),300)

q_smooth = spline(epsb,diff_5,qnew)
r_smooth = spline(epsb,diff_10,qnew)
s_smooth = spline(epsb,diff_15,qnew)
t_smooth = spline(epsb,diff_20,qnew)

'''X,Y = np.meshgrid(x,epsb)
ax.plot_wireframe(X,Y,NN_surface)
ax.plot_wireframe(X,Y,NS_surface)
ax.legend(['NN','NS'])'''
plt.plot(qnew,q_smooth,qnew,r_smooth,qnew,s_smooth,qnew,t_smooth,linewidth=3)
plt.plot(limit,np.zeros(10),'k--',linewidth=3)
plt.legend(['Radius = 5 nm','Radius = 10 nm','Radius = 15 nm','Radius = 20 nm','Zero Energy Difference'],loc=4)
plt.xlabel('Background Dielectric Constant (unitless)')
plt.ylabel('North-North/North-South Energy Difference (eV)')
plt.title('Quasistaticity as a function of background dielectric constant')
plt.xlim([1,10])
plt.show()