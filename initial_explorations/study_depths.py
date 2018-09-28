import numpy as np
import flatmaps as fm
import matplotlib.pyplot as plt
from astropy.io import fits

predir="/global/cscratch1/sd/damonge/HSC/HSC_processed/"
mask_thr=0.7
magran=[23.,27.]
snr=10
mlim=24.5
fields=['WIDE_GAMA09H','WIDE_GAMA15H','WIDE_HECTOMAP','WIDE_VVDS',
        'WIDE_WIDE12H','WIDE_AEGIS','WIDE_XMMLSS']
depth_hists={}
npixs={}
ngals={}
for field in fields :
    fsk,mp_depth=fm.read_flat_map(predir+field+'/'+field+"_%ds_depth_mean_fluxerr.npz"%snr,2)
    mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
    fskb,masked=fm.read_flat_map(predir+field+'/'+field+"_MaskedFraction.npz")
    cat=fits.open(predir+field+"/"+field+'_Catalog_i%.2lf.fits'%mlim)[1].data
    mskob=~(cat['iflags_pixel_bright_object_center']*cat['iflags_pixel_bright_object_any'])
    cat=cat[mskob]
    if fsk.get_dims()!=fskb.get_dims() :
        raise ValueError("Inconsistent fskys")
    goodpix=np.where(masked>mask_thr)[0]
    mask=np.zeros_like(masked); mask[goodpix]=1
    fsk.view_map(mp_depth*mask,posColorbar=True,
                 title=field+' $%d\\sigma$ depth'%snr,
                 xlabel='R.A.',ylabel='dec',
                 colorMin=24,colorMax=26.)
    npix=np.sum(mask*mp_depth>0)+0.
    hist,mbins=np.histogram(mask*mp_depth,bins=100,range=magran)
    depth_hists[field]=hist+0.
    npixs[field]=npix
    ngals[field]=len(cat)
    print(field,fsk.get_dims(),ngals[field],npixs[field],ngals[field]/np.sum(masked))

print(np.sum([ngals[f] for f in fields]))
plt.figure()
for field in fields :
    h=depth_hists[field]; npix=npixs[field]
    hc=np.array([np.sum(h[i:]) for i in np.arange(100)])
    plt.plot(0.5*(mbins[1:]+mbins[:-1]),hc/npix,label=field,lw=2)
h=np.sum(np.array([depth_hists[f] for f in fields]),axis=0)
npix=np.sum([npixs[f] for f in fields])
hc=np.array([np.sum(h[i:]) for i in np.arange(100)])
plt.plot(0.5*(mbins[1:]+mbins[:-1]),hc/npix,'k--',label='Total',lw=2)
plt.legend(loc='lower left')
plt.xlabel('$10\\sigma$ depth',fontsize=18)
plt.ylabel('Area fraction',fontsize=18)

pix_area=fsk.dy*fsk.dx
plt.figure()
plt.plot(0.5*(mbins[1:]+mbins[:-1]),hc*pix_area,'k-',label='Total',lw=2)
plt.xlabel('$10\\sigma$ depth',fontsize=18)
plt.ylabel('Area (deg)',fontsize=18)
plt.show()
