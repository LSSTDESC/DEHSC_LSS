import pyccl as ccl
import numpy as np
import flatmaps as fm
import matplotlib.pyplot as plt
from matplotlib import cm
from astropy.io import fits
import matplotlib

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['xtick.major.size'] = 7
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['xtick.major.pad'] = 6
matplotlib.rcParams['xtick.minor.pad'] = 6
matplotlib.rcParams['xtick.labelsize'] = 32
matplotlib.rcParams['xtick.minor.width'] = 1.6
matplotlib.rcParams['xtick.major.width'] = 1.6
matplotlib.rcParams['ytick.major.size'] = 7
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['ytick.major.pad'] = 6
matplotlib.rcParams['ytick.minor.pad'] = 6
matplotlib.rcParams['ytick.labelsize'] = 32
matplotlib.rcParams['ytick.minor.width'] = 1.6
matplotlib.rcParams['ytick.major.width'] = 1.6

times=np.array([5.705970149253732,
                6.402985074626866,
                7.629850746268656,
                8.938805970149254])
chis=3261633.44*np.array([1462.,1694.,2144.,2704.])/0.67 # lightyears
dchis=chis * np.radians(1.) * 1E-6
_,mp1=fm.read_flat_map("../data_replotting/XMMLSS/ngal_maps.fits",i_map=0)
_,mp2=fm.read_flat_map("../data_replotting/XMMLSS/ngal_maps.fits",i_map=2)
_,mp3=fm.read_flat_map("../data_replotting/XMMLSS/ngal_maps.fits",i_map=4)
_,mp4=fm.read_flat_map("../data_replotting/XMMLSS/ngal_maps.fits",i_map=6)
fsk,msk=fm.read_flat_map("../data_replotting/XMMLSS/masked_fraction.fits")
msk[msk<0.001]=0
msk[msk>=0.001]=1
cmap=cm.jet
cmap.set_bad(color='#DDDDDD')

from matplotlib.patches import Rectangle
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
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(40)
    ax.coords[1].set_major_formatter('dd')
    ax.text(0.01,0.9,r't = %.1lf billion years'%(times[i]),
            transform=ax.transAxes,fontsize=32,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='white'))
    #ax.text(0.025,0.25,r'$%d \,\,{\rm Mly}$'%(dchis[i]),
    #        transform=ax.transAxes,fontsize=32)
    #ax.add_patch(Rectangle((0.01*fsk.nx,0.07*fsk.ny),
    #                       0.15*fsk.nx,0.3*fsk.ny,
    #                       facecolor='white',alpha=0.7))
    #ax.plot([0.03*fsk.nx,0.13*fsk.nx],
    #        [0.15*fsk.ny,0.15*fsk.ny],
    #        'k-', lw=4)
    #ax.plot([0.03*fsk.nx,0.03*fsk.nx],
    #        [0.12*fsk.ny,0.18*fsk.ny],
    #        'k-', lw=4)
    #ax.plot([0.13*fsk.nx,0.13*fsk.nx],
    #        [0.12*fsk.ny,0.18*fsk.ny],
    #        'k-', lw=4)
fig.text(0.5, 0.06, 'R.A.', ha='center', fontsize=40)
fig.text(0.04, 0.5, 'Dec.', va='center', rotation='vertical', fontsize=40)
ax.coords[0].set_major_formatter('dd')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(40)
plt.savefig('../doc/galmap.pdf',bbox_inches='tight')
plt.show()
