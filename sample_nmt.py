import numpy as np
import pymaster as nmt
import flatmaps as fm
import matplotlib.pyplot as plt

DTOR=np.pi/180

#Create information for a flat-sky map covering the GAMA15H field with resolution ~1.2 arcmin
mi=fm.FlatMapInfo([212.5,222.],[-2.,2.],dx=0.024,dy=0.024)

#Generate a map as a Gaussian random field
larr=np.arange(8000.)
clarr=((larr+1000.)/1000.)**(-1.5)
mp=(nmt.synfast_flat(mi.nx,mi.ny,mi.lx*DTOR,mi.ly*DTOR,clarr)[0]).flatten()

#Generate a mask (just remove the edges of the field)
mask=np.ones([mi.ny,mi.nx])
mask[:mi.ny/5 ,:]=0; mask[ 4*mi.ny/5: ,:]=0
mask[:,:mi.nx/11]=0; mask[:,10*mi.nx/11:]=0
mask=mask.flatten()
mi.view_map(mask*mp)

#Compute power spectrum
cl,lbpw,wsp=mi.compute_power_spectrum(mp*mask,mask)

#Plot success!
plt.figure()
plt.plot(larr,clarr,'k-',label='Input $C_\\ell$')
plt.plot(np.mean(lbpw,axis=0),cl,'r-',label='Computed $C_\\ell$')
plt.xlim([lbpw[0,0],lbpw[1,-1]])
plt.xlabel('$\\ell$',fontsize=15)
plt.ylabel('$C_\\ell$',fontsize=15)
plt.legend(loc='lower left',frameon=False,fontsize=15)
plt.loglog()
plt.savefig("sample_nmt.pdf")
plt.show()
