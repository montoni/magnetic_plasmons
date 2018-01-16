import matplotlib.pyplot as plt
import numpy as np

NN_vac = np.loadtxt('NN_contour_vacuum.txt')
NS_vac = np.loadtxt('NS_contour_vacuum.txt')
NN_wat = np.loadtxt('NN_contour_water.txt')
NS_wat = np.loadtxt('NS_contour_water.txt')

vacuum = NN_vac - NS_vac
water = NN_wat - NS_wat
vaclev = np.linspace(-vacuum.max(),vacuum.max(),11) 
watlev = np.linspace(-water.max(),water.max(),11) 
plt.figure()
plt.contourf(vacuum,levels=vaclev,cmap='bwr')
plt.colorbar()
plt.title('vacuum 2mer')
plt.xlabel('separation')
plt.ylabel('scaling')
plt.savefig('vacuum_2mer_contour.pdf')
plt.show()

plt.figure()
plt.contourf(water,levels=watlev,cmap='bwr')
plt.colorbar()
plt.title('water 2mer')
plt.xlabel('separation')
plt.ylabel('scaling')
plt.savefig('water_2mer_contour.pdf')
plt.show()