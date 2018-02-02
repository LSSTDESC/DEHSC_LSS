import numpy as np
import flatmaps as fm
import matplotlib.pyplot as plt
from astropy.io import fits
from createMaps import createCountsMap

predir="/global/cscratch1/sd/damonge/HSC/HSC_processed/"
mask_thr=0.5
magran=[23.,27.]
snr=10
mlim=24.5
fields=['WIDE_GAMA09H','WIDE_GAMA15H','WIDE_HECTOMAP','WIDE_VVDS',
        'WIDE_WIDE12H','WIDE_AEGIS','WIDE_XMMLSS']
#fields=['WIDE_GAMA15H']
npix={}
cls_raw={}
cls_nd={}
cls_ns={}
cls_nds={}
for field in fields :
    print(field)
    #Create depth-based mask
    fsk,mp_depth=fm.read_flat_map(predir+field+'/'+field+"_%ds_depth_mean_fluxerr.npz"%snr,2)
    mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
    msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=mlim]=1

    #Read masked fraction
    fskb,mskfrac=fm.read_flat_map(predir+field+'/'+field+"_MaskedFraction.npz")
    if fsk.get_dims()!=fskb.get_dims() :
        raise ValueError("Inconsistent fskys")

    #Create BO-based mask
    msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>mask_thr]=1
    msk_t=msk_bo*msk_depth
    goodpix=np.where(msk_t>0.1)[0]

    #Total weight map
    weight=mskfrac*msk_t
    npix[field]=np.sum(weight)

    #Read systematics maps
    fskb,mp_dust=fm.read_flat_map(predir+field+'/'+field+"_syst_dust.npz",2)
    if fsk.get_dims()!=fskb.get_dims() :
        raise ValueError("Inconsistent fskys")
    fskb,mp_star=fm.read_flat_map(predir+field+'/'+field+"_syst_nstar_i%.2lf.npz"%mlim)
    if fsk.get_dims()!=fskb.get_dims() :
        raise ValueError("Inconsistent fskys")

    #Read catalog and cut based on BO mask
    cat=fits.open(predir+field+"/"+field+'_Catalog_i%.2lf.fits'%mlim)[1].data
    mskob=~(cat['iflags_pixel_bright_object_center']*cat['iflags_pixel_bright_object_any'])
    cat=cat[mskob]

    #Compute ndens map
    nmap=createCountsMap(cat['ra_r'],cat['dec_r'],fsk)
    ndens=np.sum(nmap*msk_t)/np.sum(mskfrac*msk_t)
    delta=np.zeros_like(weight)
    delta[goodpix]=nmap[goodpix]/(ndens*mskfrac[goodpix])-1

    #fsk.view_map(mp_depth,title='Depth '+field,colorMin=24)
    #fsk.view_map(mp_dust*weight,title='Dust '+field)
    #fsk.view_map(mp_star*weight,title='$n_{\\rm star}$ '+field)
    #fsk.view_map(delta*weight,title='$\\delta_g$ '+field)

    #Compute power spectra
    cl_raw,bpws,wsp=fsk.compute_power_spectrum(delta,weight,return_bpw=True,return_wsp=True)
    cl_nodust,bpws=fsk.compute_power_spectrum(delta,weight,return_bpw=False,return_wsp=False,
                                              wsp=wsp,temp1=[mp_dust])
    cl_nostar,bpws=fsk.compute_power_spectrum(delta,weight,return_bpw=False,return_wsp=False,
                                              wsp=wsp,temp1=[mp_star])
    cl_nodstr,bpws=fsk.compute_power_spectrum(delta,weight,return_bpw=False,return_wsp=False,
                                              wsp=wsp,temp1=[mp_dust,mp_star])
    ells=np.mean(bpws,axis=0)
    ndens_perad=ndens/(fsk.dx*fsk.dy*(np.pi/180)**2)
    cls_raw[field]=[ells,cl_raw,ndens_perad]
    cls_nd[field]=[ells,cl_nodust,ndens_perad]
    cls_ns[field]=[ells,cl_nostar,ndens_perad]
    cls_nds[field]=[ells,cl_nodstr,ndens_perad]
    #plt.figure()
    #plt.plot(ells,cl_raw-np.ones_like(ells)/ndens_perad,'r-',label='raw')
    #plt.plot(ells,cl_nodust-np.ones_like(ells)/ndens_perad,'g-',label='D')
    #plt.plot(ells,cl_nostar-np.ones_like(ells)/ndens_perad,'b-',label='S')
    #plt.plot(ells,cl_nodstr-np.ones_like(ells)/ndens_perad,'y-',label='DS')
    #plt.plot(ells,np.ones_like(ells)/ndens_perad,'k--')
    #plt.loglog()
    #plt.legend(loc='lower left')
    #plt.figure()
    #plt.plot(ells,cl_raw-cl_nodstr)
    #plt.show()
    
plt.figure()
for f in fields :
    plt.plot(cls_raw[f][0],cls_raw[f][1]-1./cls_raw[f][2],'-',label=f,lw=2)
plt.plot(cls_raw['WIDE_VVDS'][0],
         np.ones_like(cls_raw['WIDE_VVDS'][0])/cls_raw['WIDE_VVDS'][2],'k--',
         label='Shot noise')
plt.loglog()
plt.legend(loc='lower left',ncol=2,labelspacing=0.1,frameon=False)
plt.xlabel('$\\ell$',fontsize=18)
plt.ylabel('$C_\\ell$',fontsize=18)
plt.savefig('cls.pdf',bbox_inches='tight')
plt.show()
