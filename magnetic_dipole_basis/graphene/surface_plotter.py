import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
import numpy as np

x = np.linspace(-1,1,80)
[X,Y] = np.meshgrid(x,x)
Z_plus = np.loadtxt('20_nm_E_plus.txt')
Z_minus= np.loadtxt('20_nm_E_minus.txt')

ax.plot_surface(X,Y,Z_plus,color='red')
ax.plot_surface(X,Y,Z_minus,color='blue')
plt.show()
