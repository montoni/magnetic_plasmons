import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('density_0.15/epsilon_den_0.15_rad_10.0')

plt.figure()
plt.plot(data[0],data[1],data[0],data[2])
plt.legend(['real','imag'])
plt.show()