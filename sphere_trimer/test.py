import numpy as np

a = [[1 + 1j, 0.5 + 1j],[0.5 + 1j, 1 + 1j]]

w,v = np.linalg.eig(a)

print w