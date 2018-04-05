import numpy as np
import matplotlib.pyplot as plt
import math

<<<<<<< HEAD

A = np.loadtxt('../eegs_seels/eegs_point_rod_normal_2.9')[::-1]
AA = np.loadtxt('../eegs_seels/eegs_point_rod_damp_0.2_2.9')[::-1]
AAA = np.loadtxt('../eegs_seels/eegs_point_rod_damp_0.3_2.9')[::-1]
B = np.loadtxt('../eegs_seels/eels_point_rod_normal_2.9')
C = np.loadtxt('../eegs_seels/seels_point_rod_normal_2.9')
CC = np.loadtxt('../eegs_seels/seels_point_rod_damp_0.2_2.9')
CCC = np.loadtxt('../eegs_seels/seels_point_rod_damp_0.3_2.9')
D = np.loadtxt('../eegs_seels/eegs_scan_rod')
DD = np.loadtxt('../eegs_seels/ext_scan_rod')
E = np.loadtxt('../eegs_seels/eels_scan_rod')
EE = np.loadtxt('../eegs_seels/seels_scan_rod')
eV = np.linspace(1,4,151)
eV_eegs = np.flipud(eV)
plt.figure()
plt.plot(eV,D/D.max(),eV,DD/DD.max(),linewidth=3)
ax = plt.gca()
ax.invert_xaxis()
plt.legend(['rod eegs','rod extinction'])
plt.savefig('../eegs_seels/rod_scan_eegs.pdf')
plt.show()

plt.figure()
plt.plot(eV,E,eV,EE,linewidth=3)
plt.legend(['rod eels','rod seels'])
plt.savefig('../eegs_seels/rod_scan_eels.pdf')
plt.show()
'''plt.figure()
plt.plot(eV_eegs,A,eV_eegs,AA,eV_eegs,AAA,linewidth=3)
plt.legend(['eegs linewidth x1','eegs linewidth x2','eegs linewidth x3'])
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel('Energy (eV)')
plt.ylabel('Probability')
plt.savefig('../eegs_seels/rod_eegs_dark_mode_linewidth.pdf')
plt.show()
plt.figure()
plt.plot(eV,B,eV,C,eV,CC,eV,CCC,linewidth=3)
plt.legend(['eels','seels linewidth x1','seels linewidth x2','seels linewidth x3'])
plt.xlabel('Energy (eV)')
plt.ylabel('Probability')
plt.savefig('../eegs_seels/rod_eels_dark_mode_linewidth.pdf')
plt.show()'''
'''for num in range(1,31):
	print num
	A = np.loadtxt('rod_spectra/rods_'+str(num)+'nm',skiprows=1)
	B = np.loadtxt('rod_spectra/rods_NS_'+str(num)+'nm',skiprows=1)
	plt.figure()
	plt.plot(A[:,0],A[:,2])
	plt.plot(B[:,0],B[:,2])
	plt.legend(['Ey,Bz','Ey,Bx'])
	plt.show()'''

'''for num in range(1,31):
	A = np.loadtxt('reverse_bowtie/rods_NN_'+str(num)+'nm',skiprows=1)
	B = np.loadtxt('reverse_bowtie/rods_narrow_NS_'+str(num)+'nm',skiprows=1)
	plt.subplot(2,1,1)
	plt.title(str(num))
	plt.plot(A[:,0],A[:,2])
	plt.legend(['Ex,Bz'])
	plt.ylabel('Scattering Cross Section')
	plt.subplot(2,1,2)
	plt.plot(B[:,0],B[:,2])
	plt.legend(['Ey,Bx'])
	plt.xlabel('Energy (eV)')
	plt.ylabel('Scattering Cross Section')
	plt.savefig('rods_nm'+str(num)+'.pdf')
	plt.show()'''

'''for num in [1,5,10,15,20,25,30]:
	A = np.loadtxt('CL_calcs/nm'+str(num)+'_2mer_CL_updo_broad',skiprows=1)
	B = np.loadtxt('CL_calcs/nm'+str(num)+'_2mer_CL_leri_broad',skiprows=1)
	C = np.loadtxt('CL_calcs/nm'+str(num)+'_2mer_CL_foba_broad',skiprows=1)
	plt.subplot(3,1,1)
	plt.plot(A[:,0],A[:,7],label='NN up')
	plt.plot(B[:,0],B[:,7],label='NN left')
	plt.plot(C[:,0],C[:,7],label='NN for')
	plt.title(str(num)+' nm')
	plt.legend()
	plt.subplot(3,1,2)
	plt.plot(A[:,0],A[:,5],label='NS up')
	plt.plot(A[:,0],B[:,5],label='NS left')
	plt.plot(B[:,0],C[:,5],label='NS for')

	plt.subplot(3,1,3)
	plt.plot(A[:,0],A[:,9],label='x up')
	plt.plot(A[:,0],B[:,9],label='x left')
	plt.plot(B[:,0],C[:,9],label='x for')
	#plt.title(str(num)+' nm')
	plt.legend()
	plt.show()'''

'''for num in [1,15,30,40,50]:
	A = np.loadtxt('CL_equal_hemi/nm'+str(num)+'_2mer_CL_updo',skiprows=1)
	B = np.loadtxt('CL_equal_hemi/nm'+str(num)+'_2mer_CL_leri',skiprows=1)
	C = np.loadtxt('CL_equal_hemi/nm'+str(num)+'_2mer_CL_foba',skiprows=1)
	plt.subplot(3,1,1)
	plt.plot(A[:,0],A[:,5],label='NS up')
	plt.plot(A[:,0],A[:,6],label='NS down')
	#plt.plot(C[:,0],C[:,5],label='NS for')
	#plt.plot(B[:,0],B[:,5]/(B[:,5].max()),label='NS left')
	plt.title(str(num)+' nm')
	plt.legend()
	plt.subplot(3,1,2)
	plt.plot(A[:,0],B[:,7],label='NN left')
	plt.plot(A[:,0],B[:,8],label='NN right')
	#plt.plot(C[:,0],C[:,7],label='NN for')
	
	plt.legend()
	plt.subplot(3,1,3)
	plt.plot(A[:,0],C[:,7],label='NN for')
	plt.plot(A[:,0],C[:,8],label='NN back')
	#plt.plot(C[:,0],C[:,9],label='both for')
	plt.legend()
	plt.show()'''
'''phi = np.linspace(0,2*math.pi,181)
for num in [20]:
	for ene in np.linspace(3.3,3.5,21):
		if ene == 3.0:
			ene = int(ene)
		elif ene == 4.0:
			ene = int(ene)
		else:
			pass
		A = np.loadtxt('CL_angular/yz_NS_CL_'+str(num)+'_ene_'+str(ene))
		B = np.loadtxt('CL_angular/xz_NN_CL_'+str(num)+'_ene_'+str(ene))
		#f, axarr = plt.subplots(1,1,subplot_kw=dict(projection='polar'))
		plt.figure()
		
		plt.polar(phi,A,label='NS yz')
		plt.polar(phi,B,label='NN xz')
		#axarr[1] = plt.plot(phi,B,label='xz')
		plt.title(str(num) + 'nm, ' +str(ene)+ ' eV')
		plt.legend()
		plt.show()'''
=======
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
>>>>>>> 248f1ce6598f71d449290d94abddb91ea29a5b3d
