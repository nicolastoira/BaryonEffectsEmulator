#Define path to BaryonEffectsEmulator main code and import its functions
import sys
sys.path.insert(0, '/path/to/BaryonEffectsEmulator-master/Main source')
from BaryonEffectsEmulator import *

#Define baryonic collection model
MyBCM = {'f_b': 0.14, 'logMc': 14, 'mu': 0.4, 'theta_ej':4.0, 'eta_tot': 0.3, 'eta_cga': 0.6}

#Define input redshifts
z=[0.3,0.9,1.5]

#Define wavenumbers of evaluation
keval=np.logspace(-1,1,1000)

#Calculate the boosts
result=get_boost(z,MyBCM,keval)

#Print some results
print 'k-values: \n'+str(result['k'])+'\n'
print 'Boost for redshift z= '+str(z[0])+':\n'+str(result['z0'])

#Plot the ratios
plt.figure(figsize=(10,5))
plt.semilogx(result['k'],result['z0'],label=r'$z=$'+str(z[0]))
plt.semilogx(result['k'],result['z1'],label=r'$z=$'+str(z[1]))
plt.semilogx(result['k'],result['z2'],label=r'$z=$'+str(z[2]))
plt.axhline(y=1,color='k')
plt.title(r'Emulation of the power spectrum ratio for a given redshift and BCM')
plt.xlabel(r'$k$ $[h/\mathrm{Mpc}]$')
plt.ylabel(r'P$_{\rm{BCM}}$ /P$_{\rm{DMO}}$')
plt.legend(loc='best',frameon=False)
plt.grid()
plt.savefig('example_plot.png',dpi=500)