import numpy as np
import flatmaps as fm
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits


_,mp1=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/ngal_maps.fits",i_map=0)
_,mp2=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/ngal_maps.fits",i_map=2)
_,mp3=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/ngal_maps.fits",i_map=4)
_,mp4=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/ngal_maps.fits",i_map=6)
fsk,msk=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/masked_fraction.fits")
msk[msk<0.001]=0
msk[msk>=0.001]=1
#fig=plt.figure(figsize=(20,fsk.ny*20./fsk.nx))
#ax=fig.add_subplot(111,projection=fsk.wcs)
#im=ax.imshow(msk.reshape([fsk.ny,fsk.nx]),vmin=0,vmax=1,
#             cmap=cm.binary, origin='lower', interpolation='nearest')
cmap=cm.jet
cmap.set_bad(color='#AAAAAA')

from scipy.ndimage import gaussian_filter
fig=plt.figure(figsize=(20,4*fsk.ny*20./fsk.nx))
plt.subplots_adjust(hspace=0)
for i,mp in enumerate([mp1,mp2,mp3,mp4]):
    ax=fig.add_subplot(4,1,i+1,projection=fsk.wcs)
    mean=np.sum(mp*msk)/np.sum(msk)
    delta=msk*(mp-mean)/mean
    delta=gaussian_filter(delta.reshape([fsk.ny,fsk.nx]),sigma=1).flatten()
    delta[msk<1]=np.nan
    im=ax.imshow(delta.reshape([fsk.ny,fsk.nx]),vmin=-1,vmax=4.5,
                 cmap=cmap, origin='lower', interpolation='nearest')
    ax.set_ylabel('${\\rm Dec.}\\,(^\\circ)$',
                  fontsize=14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
ax.set_xlabel('${\\rm R.A.}\\,(^\\circ)$',
              fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
#    fig.colorbar(im)
plt.savefig('galmap.pdf',bbox_inches='tight')
plt.show()
exit(1)
cat=fits.open("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/clean_catalog.fits")[1].data
z0=0.
zf=10
cat=cat[cat['pz_best_eab']<zf]
cat=cat[cat['pz_best_eab']>z0]
ix,iy=np.transpose(fsk.wcs.wcs_world2pix(np.transpose(np.array([cat['ra'],
                                                                cat['dec']])),0))
ax.scatter(ix,iy,c=(cat['pz_best_eab']-z0)/(zf-z0),s=1,cmap=cm.viridis,edgecolors='none')
print(len(cat))
print(cat.names)
plt.show()
