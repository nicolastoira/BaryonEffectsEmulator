
import numpy as np
import pandas as pd
from scipy import special
import os
import matplotlib.pyplot as plt
import time
import math

"""
Function get_boost(z,BCM_params) returns the k-values and the evaluated boost factor for a chosen redshift z=0,0.5,1,1.5,2

INPUT: redshift z, BCM free parameters

OUTPUT: k-values, boost factor for each k-value

"""

def get_boost(z,BCM_params):
    nSteps = 4
    nk = 232
    nz = nSteps+1
    #omega_rad = 4.183709424969527e-5; /* corresponds to 2.755 K Tcmb */
    #bVerbose = 0;
    nCoef = np.array([1340,1706,1706,1706,1706,1706,1706,1706,1340,1706,1706,1706,1340,1340,1340,1706,1706,1706,1706,1340,1064,219,1706,456])

    min_par = np.array([0.13,12.7,0.1,2.0,0.2,0.5])
    max_par = np.array([0.21,16.7,1.0,8.0,0.4,0.7])
    nPCA=24
    nPar=6

    start=time.time()

    #fd = pd.read_csv("./ee_bar_new.dat", sep=' ',header=None)
    #fd.to_pickle("./ee_bar_new.pkl")
    fd=pd.read_pickle("/Users/Nicola/Desktop/Thesis/Data/wrapper/ee_bar_new.pkl")
    f_val=fd.values

    
    # Reading in PCE coefficients
    coef=np.zeros(np.sum(nCoef))
    di=0
    for i in range(np.sum(nCoef)):
        coef[i]=f_val[0][i]
    di=di+np.sum(nCoef)

    #Reading in PCE multi-indices
    index=np.zeros(np.sum(nPar*nCoef))
    for i in range(np.sum(nPar*nCoef)):
        index[i]=f_val[0][di+i]
    di=di+np.sum(nPar*nCoef)


    #Reading in principal components
    pca=np.zeros(nPCA*nk*nz)
    for i in range(nPCA*nk*nz):
        pca[i]=f_val[0][di+i]
    di=di+nPCA*nk*nz


    #Reading in PCA means
    mean=np.zeros(nk*nz)
    for i in range(nk*nz):
        mean[i]=f_val[0][di+i]
    di=di+nk*nz


    #Reading in k vector 
    kvec=np.zeros(nk)
    for i in range(nk):
        kvec[i]=f_val[0][di+i]
    di=di+nk

    #Determine lmax (= maximal order of Legendre polynomials)
    lmax=0
    for i in range(len(index)):
        l=index[i]
        if l>lmax:
            lmax=l
    
    #Calculate all the needed Legendre polynomials up to order lmax for each of the cosmological input parameters (scaled to [-1,1)).
    Pl=np.zeros((6,int((lmax+1))))

    for ip in range(6):
        x=2*(BCM_params[ip]-min_par[ip])/(max_par[ip]-min_par[ip])-1.0
        Pl[ip][:]=special.eval_legendre(np.arange(0,lmax+1), x)
        for ll in range(int(lmax)+1):
            Pl[ip][ll]=Pl[ip][ll]*np.sqrt(2.0*ll + 1.0)

    #Now assemble the eigenvalues for the PCAs
    lamb=np.zeros(24)
    number=np.array([0,1,2,3,4,5])
    for i in range(24):
        for ic in range(nCoef[i]):
            idx=index[6*np.sum(nCoef[0:i])+number[:]*nCoef[i]+ic].astype(int)
            prod=np.prod(Pl[number[:],idx[number[:]]])
            prod=prod*coef[np.sum(nCoef[0:i])+ic]
            lamb[i]=lamb[i]+prod
  
    """ Original version a bit slower than the new one
    for i in range(24):
        lamb[i]=0
        for ic in range(nCoef[i]):
            prod=1.0
            for ip in range(6):
                prod=prod*Pl[ip][int(index[6*np.sum(nCoef[0:i])+ip*nCoef[i]+ic])]
            prod=prod*coef[np.sum(nCoef[0:i])+ic]
            lamb[i]=lamb[i]+prod
    """       

    #Boost
    len_z=1 #provvisory
    boost=np.zeros(len_z*nk)
    lnboost=np.zeros(len_z*nk)
    zcounter=0

    #Following to lines to uncomment if interpolation works
    #if z==0:
    #   iz=0
    zs=[0.,0.5,1.,1.5,2]

    if z in zs:
        for i in range(nz):
            if z==zs[i]:
                iz=i

        for j in range(nk):
            lnboost[zcounter*nk+j]=mean[iz*nk+j]

        for j in range(nk):
            for i in range(24):
                lnboost[zcounter*nk+j]=lnboost[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+iz*nk+j]


    #Interpolation for an arbitrary redshift    
    
    
    else:
        boost0=np.zeros(len_z*nk)
        lnboost0=np.zeros(len_z*nk)
        boost1=np.zeros(len_z*nk)
        lnboost1=np.zeros(len_z*nk)

        min_dist=min(zs, key=lambda variable:abs(variable-z))
        if min_dist-z>0:
            z0=min_dist-0.5
            z1=min_dist
        else:
            z0=min_dist
            z1=min_dist+0.5

        for i in range(nz):
            if z0==zs[i]:
                iz0=i
            if z1==zs[i]:
                iz1=i

        for j in range(nk):
            lnboost0[zcounter*nk+j]=mean[iz0*nk+j]
            lnboost1[zcounter*nk+j]=mean[iz1*nk+j]

        for j in range(nk):
            for i in range(24):
                lnboost0[zcounter*nk+j]=lnboost0[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+iz0*nk+j]
                lnboost1[zcounter*nk+j]=lnboost1[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+iz1*nk+j]

        lnboost=lnboost0+(lnboost1-lnboost0)/(z1-z0)*(z-z0)

    #USING SPLINE
    """
    else:
        boost0=np.zeros(len_z*nk)
        lnboost0=np.zeros(len_z*nk)
        boost05=np.zeros(len_z*nk)
        lnboost05=np.zeros(len_z*nk)
        boost1=np.zeros(len_z*nk)
        lnboost1=np.zeros(len_z*nk)
        boost15=np.zeros(len_z*nk)
        lnboost15=np.zeros(len_z*nk)
        boost2=np.zeros(len_z*nk)
        lnboost2=np.zeros(len_z*nk)

        for j in range(nk):
            lnboost0[zcounter*nk+j]=mean[0*nk+j]
            lnboost05[zcounter*nk+j]=mean[1*nk+j]
            lnboost1[zcounter*nk+j]=mean[2*nk+j]
            lnboost15[zcounter*nk+j]=mean[3*nk+j]
            lnboost2[zcounter*nk+j]=mean[4*nk+j]

        for j in range(nk):
            for i in range(24):
                lnboost0[zcounter*nk+j]=lnboost0[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+0*nk+j]
                lnboost05[zcounter*nk+j]=lnboost05[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+1*nk+j]
                lnboost1[zcounter*nk+j]=lnboost1[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+2*nk+j]
                lnboost15[zcounter*nk+j]=lnboost15[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+3*nk+j]
                lnboost2[zcounter*nk+j]=lnboost2[zcounter*nk+j]+lamb[i]*pca[nk*nz*i+4*nk+j]
        from scipy.interpolate import splev, splrep
        for i in range(nk*len_z):
            spl=splrep(zs,[lnboost0[i],lnboost05[i],lnboost1[i],lnboost15[i],lnboost2[i]])
            lnboost[i]=splev(z,spl)
    """

    for j in range(len_z*nk):
            boost[j]=np.exp(lnboost[j])
    end=time.time()
    print 'Runtime: ' +str(end-start)+' s'
    return kvec,boost
    


"""
Fucntion check_param_range(par_dict) checks if the input parameters obey the limits given by the emulator

INPUT: dictionary containing the free parameters values

OUTPUT: errors if one or more of the parameters do not obey the limits of the emulator

"""

def check_param_range(par_dict,z,csm_index=0):
    """
    Checks if all parameters in the cosmology dictionary 'par_dict'
    (with index 'csm_index') passed to this function obey the limits
    set by the emulator.

    """
    f_b = [0.13, 0.21]
    logMc = [12.7, 16.7]
    mu = [0.1, 1.0]
    theta_ej = [2.0,8.0]
    eta_tot = [0.2, 0.4]
    eta_cga = [0.5, 0.7]
    z_values=[0.,2.]

    f_b_not_in_range = f_b[0] > par_dict['f_b'] or f_b[1] < par_dict['f_b']

    logMc_not_in_range = logMc[0] > par_dict['logMc'] or logMc[1] < par_dict['logMc']

    mu_not_in_range = mu[0] > par_dict['mu'] or mu[1] < par_dict['mu']

    theta_ej_not_in_range = theta_ej[0] > par_dict['theta_ej'] or theta_ej[1] < par_dict['theta_ej']

    eta_tot_not_in_range = eta_tot[0] > par_dict['eta_tot'] or eta_tot[1] < par_dict['eta_tot']

    eta_cga_not_in_range = eta_cga[0] > par_dict['eta_cga'] or eta_cga[1] < par_dict['eta_cga']

    z_not_in_range = z_values[0] > z or z_values[1] < z 
   

    if f_b_not_in_range:
        raise ValueError("Parameter range violation: f_b is set to %f, \
                          but should be in the interval [0.13, 0.21]." \
                      %(par_dict['f_b']))

    if logMc_not_in_range:
        raise ValueError("Parameter range violation: logMc is set to %f, \
                          but should be in the interval [12.7, 16.7]." \
                         %(par_dict['logMc']))

    if mu_not_in_range:
        raise ValueError("Parameter range violation: mu is set to %f, \
                          but should be in the interval [0.1, 1.0]." \
                         %(par_dict['mu']))

    if theta_ej_not_in_range:
        raise ValueError("Parameter range violation: theta_ej is set to %f, \
                          but should be in the interval [2.0,8.0]." \
                         %(par_dict['theta_ej']))

    if eta_tot_not_in_range:
        raise ValueError("Parameter range violation: eta_tot is set to %f, \
                          but should be in the interval [0.2, 0.4]." \
                         %(par_dict['eta_tot']))

    if eta_cga_not_in_range:
        raise ValueError("Parameter range violation: eta_cga is set to %f, \
                          but should be in the interval [0.5, 0.7]." \
                         %(par_dict['eta_cga']))

    if z_not_in_range:
        raise ValueError("Parameter range violation: z is set to %f, \
                          but should be in the interval [0.0,2.0]." \
                         %z)


def BarEmu_plot(dictionary,z):
    #z=raw_input('Redshift: ')

    #Check parameters
    check_param_range(dictionary,z)

    BCM_parameters=[dictionary['f_b'],dictionary['logMc'],dictionary['mu'],dictionary['theta_ej'],dictionary['eta_tot'],dictionary['eta_cga']]
   
    #Emulate the boost for the given BCM parameters
    kvec,boost=get_boost(float(z),BCM_parameters)


    fig=plt.figure(figsize=(8,5))
    plt.semilogx(kvec,boost,'r-',label=r'Emulated boost at $z=$'+str(z))
    plt.annotate(BCM_parameters, xy=(0.1, 0.1), xycoords='axes fraction')
    plt.axhline(y=1, color='k', linestyle='-')
    plt.xlabel(r'$k$ $[h/\mathrm{Mpc}]$')
    plt.ylabel(r'P$_{\rm{BCM}}$ /P$_{\rm{DMO}}$')
    plt.legend(loc='best',frameon=False)
    plt.xlim(1e-1,20)
    plt.grid()
    plt.title(r'P$_{\rm{BCM}}$ /P$_{\rm{DMO}}$ logarithmic emulation for input model')
    #plt.show()
    return kvec,boost,fig

