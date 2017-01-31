import numpy as np
import compute
import math
#import os
#os.system("rm -r _output")
#os.system("mkdir _output")


# Choose dot radii and set density values
R = np.linspace(1.5,30.,58)
#R = 4.0
den = 0.15 # nm^-3
print R
# Compute nF's from shell model (could make this another script
# that ouputs a text file?)
vol = (4./3.)*(np.pi)*(R**3)
M = np.floor(0.5*vol*den)
print M
nF = []
for m in M:
    coeff = [1./3., 1./2., 1./6., -1.*m]
    a = np.roots(coeff)
    #print a
    nF.append(int(np.around(a[-1].real)))
print nF

# Create list of names
names = []
for r in R:
    names.append('den_' + str(den) + '_rad_' + str(r))


for i in xrange(0,len(R)):
    compute.eps(R[i],nF[i],den,names[i])

