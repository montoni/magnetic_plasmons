p_dot_e = 0
for x in range(0,numPart):
    electric = np.multiply([0,1],np.exp(1j*wavenumber*np.sqrt(Loc[x][0]**2 + Loc[x][1]**2)))
p_dot_E += np.dot(vec[x],electric)
screen_dist = 1000000000000*a0
vec = vec*p_dot_E
for number in range(0,numPart):
    location = list(Loc[number])
location.append(0)

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
Bfield[theta_count,phi_count] = (wavenumber**2)*np.cross(nhat_dip,vec[number])*np.exp(1j*wavenumber*rmag)/rmag #
Efield[theta_count,phi_count] = np.cross(Bfield[theta_count,phi_count],nhat_dip)/(epsb)
theta_count += 1
phi_count += 1
Bfield_total = Bfield_total + Bfield
Efield_total = Efield_total + Efield
poynting = np.empty((poynting_points,poynting_points),dtype=float)
for idx1 in range(poynting_points):
    for idx2 in range(poynting_points):
        poynting[idx1,idx2] = (screen_dist**2)*np.linalg.norm(np.real(np.cross(Efield_total[idx1,idx2],np.conj(Bfield_total[idx1,idx2]))))
theta = np.linspace(0,math.pi,poynting_points)
phi = np.linspace(0,2*math.pi,poynting_points)
PHI,THETA = np.meshgrid(phi,theta)
X = poynting * np.cos(PHI)*np.sin(THETA)
Y = poynting * np.sin(PHI)*np.sin(THETA)
Z = poynting * np.cos(THETA)
norm = poynting/poynting.max()
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm),linewidth=0, antialiased=False, alpha=0.5)
plt.show()
