import numpy as np
import matplotlib.pyplot as plt
import flatmaps as fm
from scipy.interpolate import interp1d
import formatting

# Read best-fit Cls and interpolate them
d=np.loadtxt("../hsc_lss_params/cls_guess_i24p5.txt",unpack=True)
# Read assuming constant Cls for l<2
ls=np.arange(int(d[0][-1])+1,dtype=float)
cls=np.zeros([4,len(ls)])
for i,i_d in enumerate([1,5,8,10]):
    cls[i,2:]=d[i_d]
cls[:,:100]=cls[:,100][:,None]
lclf=[interp1d(ls,np.log10(c),bounds_error=False,fill_value=np.log10(c[-1])) for c in cls]

def get_ndens_and_error(field_name):
    # This computes the number density and associated uncertainty
    # for all redshift bins for a given field.

    # Read mask
    prefix="/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+field_name+"_sirius_i24p5_out/"
    prefix="../data_replotting/"+field_name+"/"
    fsk,mskfrac=fm.read_flat_map(prefix+"masked_fraction.fits")
    msk_bin=np.zeros_like(mskfrac)
    msk_bin[mskfrac>0.5]=1
    mskfrac*=msk_bin;

    # Footprint area (in arcmin)
    area=np.sum(mskfrac)*fsk.dx*fsk.dy*60**2
    
    # Count number of galaxies in each redshift bin
    ngals=np.array([np.sum(msk_bin*
                           fm.read_flat_map(prefix+"ngal_maps.fits",
                                            i_map=2*i)[1])
                    for i in range(4)])

    # Ell sampling for DFTs
    dlx=2*np.pi/np.radians(fsk.lx)
    dly=2*np.pi/np.radians(fsk.ly)
    ell_y=np.fft.fftfreq(fsk.ny,d=1./(dly*fsk.ny))
    ell_x=np.fft.fftfreq(fsk.nx,d=1./(dlx*fsk.nx))
    # |l|
    ell_mod=np.sqrt(ell_y[:,None]**2+ell_x[None,:]**2)
    # Mask Fourier transform
    w2_fourier=np.absolute(np.fft.fft2(mskfrac.reshape([fsk.ny,fsk.nx]))/np.sum(mskfrac))**2
    # Var(delta) = Int[dl^2 * |W(l)|^2 * C_l]/(2*pi)^2
    sigma_delta=np.array([np.sqrt(dlx*dly/(2*np.pi)**2 * #dl_y * dl_x/2pi^2
                                  np.sum(w2_fourier*10.**lc(ell_mod))) # Sum[W^2 * C_ell]
                          for lc in lclf])

    # Number densities
    ndens_amin=ngals/area
    # Poisson error
    sigma_ndens_sn=np.sqrt(ngals)/area
    # Total error
    sigma_ndens_cv=sigma_delta*ngals/area

    return area, ndens_amin, sigma_ndens_cv, sigma_ndens_sn

# Compute densities and uncertainties for each field and bin
field_names=['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']
nfields=len(field_names)
ifields=np.arange(nfields)
areas=[]
ndenss=[]
sigcvs=[]
sigsns=[]
for f in field_names:
    a,n,sc,ss=get_ndens_and_error(f)
    areas.append(a)
    ndenss.append(n)
    sigcvs.append(sc)
    sigsns.append(ss)
areas=np.array(areas)
ndens=np.array(ndenss).T
sigcvs=np.array(sigcvs).T
sigsns=np.array(sigsns).T
ndens_mean=np.sum(ndens/sigcvs**2,axis=1)/np.sum(1./sigcvs**2,axis=1)
print(np.amax(np.fabs((ndens-ndens_mean[:,None])/sigcvs)))

# Plot results
fig,axes=plt.subplots(4,1,figsize=(8,8),sharex=True)
plt.subplots_adjust(left=0.2)
for ib,(ax,n,nm,s) in enumerate(zip(axes,ndens,ndens_mean,sigcvs)):
    ax.errorbar(ifields,n,yerr=s,fmt='.',c=formatting.cmap1[2],markersize=13, elinewidth=2.2, capthick=2.2, capsize=3.5,)
    ax.plot([-1,nfields],[nm,nm],'k--')
    ax.set_xticks(ifields)
    ax.set_xticklabels([])
    ax.tick_params(labelsize="large")
    ax.text(0.9,0.8,'Bin %d'%(ib+1),transform=ax.transAxes,fontsize=18)
    ax.set_xlim([-0.5,nfields-0.5])
    ax.set_ylabel('$\\bar{n}\\,[{\\rm arcmin}^{-2}]$',fontsize=18)#,va='center',rotation='vertical')
ax.set_xticklabels(field_names)
plt.savefig('../doc/Paper/figures/ndens_consistency.pdf',
            bbox_inches='tight')
plt.show()
