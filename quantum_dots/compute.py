def eps(R,nF,den,name):
    import math
    import numpy as np
    # physical parameters
    M=300
    pi=math.pi
    w_P = 0.928317 #[eV/hbar]
    hbar = 6.582119e-16 #[eVs]
    hbarst = 1.05457e-34 #[Js]
    Me = 0.24*9.10938e-31 #[kg] electron mass

    # Set the rest of the physical parameters
    Vol = (4.0/3.0)*pi*(R)**3 #[nm^3] QD volume
    gamma = 0.075+(0.1*0.553)/R #[eV/hbar] gamma
    coeff = (0.391697)/(R**2)
    
    # Loop to initial omegas
    omegas = np.linspace(0.1,6,591) # these are in eV
    epsinf = 3.71 #epsilon infinity for ZnO
    epsilon_real = np.zeros(np.size(omegas))
    epsilon_imag = np.zeros(np.size(omegas))
    for i in range(0,np.size(omegas)):

        w = omegas[i]
      
        for dl in [-1,1]:
            for dn in range(((1-dl)/2),nF+1):
                for n in range(0,(1+nF-(1-dl)/2)):
                    
                    low = 2*(nF-n-dn+(1-dl)/2)
                    
                    for l in range(low,(2*(nF-n)+1)):
                    
                        S = S_if(nF,n,dn,l,dl)
                        w_if = coeff*(4.0*n+2.0*dn+2.0*l+dl+4.0)*(2.0*dn+dl)
                        denom = (w_if**2 - w**2)**2 + (gamma*w)**2
                        epsilon_real[i] = epsilon_real[i] + (S*(w_P**2)*(w_if**2 - w**2))/denom  #real
                        epsilon_imag[i] = epsilon_imag[i] + (S*gamma*w*(w_P**2))/denom
                    
        
    epsilon_real = epsilon_real + epsinf 
    #print epsilon_real
    #print epsilon_imag
    #epsilon = [epsilon_real, epsilon_imag]
    
    data_to_save = [omegas, epsilon_real, epsilon_imag]
    np.savetxt('_'.join(['epsilon',str(name)]),data_to_save)

    return

def S_if(nF,n,dn,l,dl):
    import math
    import numpy as np
    pi = math.pi
    
    term1 = 0.0
    term2 = 0.0
    
    if dl == -1:
        denom = ((pi**2)*(nF**3)*((4.0*n+2.0*dn+2.0*l+3.0)**3)*((2.0*dn-1.0)**3))
        
        term1 = (16.0*l*((2.0*n+2.0*dn+l+1)**2)*((2.0*n+l+2.0)**2))/denom
    else:
        pass
    
    
    if dl == 1:
        denom = ((pi**2)*(nF**3)*((4.0*n+2.0*dn+2.0*l+5.0)**3)*((2.0*dn+1.0)**3))

        term2 = (16.0*(l+1)*((2.0*n+2.0*dn+l+3.0)**2)*((2.0*n+l+2.0)**2))/denom
    else: 
        pass
    
    S_if = term1 + term2
    return S_if
