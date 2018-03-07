import numpy as np
import matplotlib.pyplot as plt

'''for rad in [5]:
	for eV in np.linspace(3,3.8,41):
		if eV == 3.0:
			eV = 3
		else:
			pass
		NN_field = np.loadtxt('fields/2mer_nn/NN_2merBfield_'+str(rad)+'_nm_'+str(eV))
		NS_field = np.loadtxt('fields/2mer_ns/NS_2merBfield_'+str(rad)+'_nm_'+str(eV))
		NN_field = np.nan_to_num(NN_field)
		NS_field = np.nan_to_num(NS_field)
		Nlevels = np.linspace(NN_field.min(),NN_field.max(),21)
		Slevels = np.linspace(NS_field.min(),NS_field.max(),21)
		plt.subplot(2,1,1)
		plt.contourf(NN_field,levels=Nlevels,cmap='bwr')
		plt.colorbar()
		plt.title('NN '+str(eV))
		plt.subplot(2,1,2)
		plt.contourf(NS_field,levels=Slevels,cmap='bwr')
		plt.colorbar()
		plt.title('NS '+str(eV))
		plt.show()'''

'''80:-1'''
'''0:101'''
'''eV = np.linspace(3.3,3.6,31)'''

'''A = np.loadtxt('../eegs_seels/eegs_rod_3.3')
AA = np.loadtxt('../eegs_seels/eels_rod_3.3')
AAA = np.loadtxt('../eegs_seels/seels_rod_3.3')

eV = np.linspace(1,4,151)

plt.figure()
plt.plot(eV,AA,eV,AAA)
plt.legend(['EELS','SEELS'])
plt.xlabel('Energy (eV)')
plt.ylabel('Probability')
plt.savefig('../eegs_seels/eels_seels_rod_3.3.pdf')
plt.show()'''

'''n = 1
eV = 1.6
levels = np.linspace(0,1,21)
A = np.loadtxt('../eegs_seels/n='+str(n)+'/rod_EEGS_map_'+str(eV))
B = np.loadtxt('../eegs_seels/n='+str(n)+'/rod_EELS_map_'+str(eV))
C = np.loadtxt('../eegs_seels/n='+str(n)+'/rod_SEELS_map_'+str(eV))
plt.figure()
plt.contourf(A,levels=np.linspace(A.min(),A.max(),21))
plt.colorbar()
plt.savefig('../eegs_seels/rod_eegs_n='+str(n)+'_'+str(eV)+'.pdf')
plt.show()
EELS = np.concatenate([B[0:len(B)/2-1],B[len(B)/2+1:-1]],axis=1)
SEELS = np.concatenate([C[0:len(C)/2-1],C[len(C)/2+1:-1]],axis=1)
diff = SEELS-EELS
normdiff = diff/diff.max()
plt.figure()
plt.contourf(EELS,levels=np.linspace(EELS.min(),EELS.max(),21))
plt.colorbar()
plt.savefig('../eegs_seels/rod_eels_n='+str(n)+'_'+str(eV)+'.pdf')
plt.show()

plt.figure()
plt.contourf(SEELS,levels=np.linspace(SEELS.min(),SEELS.max(),21))
plt.colorbar()
plt.savefig('../eegs_seels/rod_seels_n='+str(n)+'_'+str(eV)+'.pdf')
plt.show()

plt.figure()
plt.contourf(diff,levels=np.linspace(diff.min(),diff.max(),21))
plt.colorbar()
plt.savefig('../eegs_seels/rod_diff_n='+str(n)+'_'+str(eV)+'.pdf')
plt.show()'''


for num in [50]:
	#B = np.loadtxt('CL_points/nm'+str(num)+'_anth_CL_NSN',skiprows=1)
	A = np.loadtxt('greybush_2017/nm'+str(num)+'_13_CL_NN',skiprows=1)
	C = np.loadtxt('greybush_2017/nm'+str(num)+'_13_CL_NS',skiprows=1)
	#plt.subplot(2,1,1)
	plt.figure()
	plt.plot(A[:,0],A[:,2],C[:,0],C[:,4])
	plt.legend(['NN xz','NS z'])
	#plt.subplot(2,1,2)
	#plt.plot(A[80:-1,0],B[80:-1,1],A[80:-1,0],B[80:-1,4])
	#plt.legend(['NN xz', 'NN yz'])
	#plt.subplot(3,1,3)
	#plt.plot(C[:,0],C[:,1],C[:,0],C[:,4])
	#plt.legend(['x xz', 'x yz'])
	#plt.savefig('greybush_2017/31_CL_'+str(num)+'nm.pdf')
	plt.show()

'''for num in [1]:
	A = np.loadtxt('CL_points/nm15_phen_CL_NN1',skiprows=1)
	B = np.loadtxt('CL_points/nm15_phen_CL_NS2',skiprows=1)
	C = np.loadtxt('CL_points/nm30_phen_CL_NN1',skiprows=1)
	D = np.loadtxt('CL_points/nm30_phen_CL_NS2',skiprows=1)
	plt.subplot(2,1,1)
	#plt.figure()
	plt.plot(A[80:-1,0],A[80:-1,2]/A[80:-1,2].max(),B[80:-1,0],B[80:-1,2]/B[80:-1,2].max())
	plt.legend(['NN yz','NS xz'])
	plt.subplot(2,1,2)
	plt.plot(C[0:121,0],C[0:121,2]/C[0:121,2].max(),D[0:121,0],D[0:121,2]/D[0:121,2].max())
	plt.legend(['NN yz', 'NS xz'])
	#plt.subplot(3,1,3)
	#plt.plot(C[:,0],C[:,1],C[:,0],C[:,4])
	#plt.legend(['x xz', 'x yz'])
	plt.savefig('phen_CL_both.pdf')
	plt.show()'''

'''for num in range(0,6):
	A = np.loadtxt('bowtie/rods_spacing_NS_sep'+str(num),skiprows=1)
	B = np.loadtxt('bowtie/rods_spacing_NN_sep'+str(num),skiprows=1)
	plt.subplot(2,1,1)
	plt.plot(A[:,0],A[:,2])
	plt.subplot(2,1,2)
	plt.plot(B[:,0],B[:,2])
	plt.legend(['Ey,Bx','Ex,Bz'])
	#plt.savefig('3mer_spectra/vac_phen_upright_'+str(num)+'nm.pdf')
	plt.show()'''

'''for num in [1]:
	A = np.loadtxt('greybush_2017/30nm_thirtyone_greybush',skiprows=1)
	plt.plot(A[:,0],A[:,2],label='mag abs')
	plt.plot(A[:,0],A[:,4],label='ele abs')
	#plt.plot(A[:,0],A[:,3],label='ele ext')
	#plt.plot(A[:,0],A[:,4],label='ele sca')
	plt.legend()
	plt.show()'''

'''for num in [1,5,10,15,20,25,30]:
	A = np.loadtxt('point_spectra/yx_vac_twomer_points_sca_'+str(num),skiprows=1)
	B = np.loadtxt('point_spectra/yz_vac_twomer_points_sca_'+str(num),skiprows=1)
	plt.plot(1240./A[:,0],A[:,1],label='mag for')
	plt.plot(1240./B[:,0],B[:,1],label='ele up')
	plt.legend()
	plt.show()'''

'''for num in [10,20,30,40,50]:
	A = np.loadtxt('greybush_CL/au_nm'+str(num)+'_31_CL_NN',skiprows=1)
	B = np.loadtxt('greybush_CL/au_nm'+str(num)+'_31_CL_NS',skiprows=1)
	#plt.subplot(2,1,1)
	plt.figure()
	plt.plot(A[:,0],A[:,1]/A[:,1].max(),B[:,0],B[:,1]/B[:,1].max())
	plt.legend(['NN yz','NS xz'])
	#plt.subplot(2,1,2)
	#plt.subplot(3,1,3)
	#plt.plot(C[:,0],C[:,1],C[:,0],C[:,4])
	#plt.legend(['x xz', 'x yz'])
	plt.savefig('phen_CL_both.pdf')
	plt.show()'''