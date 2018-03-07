import numpy as np
import matplotlib.pyplot as plt
import math


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