import flatmaps as fm
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

prefix="/global/cscratch1/sd/damonge/HSC_ceci/WIDE_XMMLSS_sirius_i24p5_out/"
fsk,msk=fm.read_flat_map(prefix+"masked_fraction.fits")
msk_bin=np.zeros_like(msk); msk_bin[msk>0.]=1

def plot_syst(name,fname,band,units=None,savename=None,scale=1):
    label=name
    if band is None:
        imap=0
    else:
        imap=['g','r','i','z','y'].index(band)
        label+=" "+band
    if units is not None:
        label+=" ["+units+"]"

    _,mp=fm.read_flat_map(prefix+fname,i_map=imap)
    mp_mean=np.sum(msk_bin*msk*mp)/np.sum(msk_bin*msk)
    mp_plot=(msk_bin*(mp-mp_mean)+mp_mean)*scale
    mp_min=np.amin(mp_plot)
    mp_max=np.amax(mp_plot)
    if name == "Star count":
        mp_max=4
    mp_plot[msk_bin<1]=mp_min-1
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=fsk.wcs)
    im=ax.imshow(mp_plot.reshape([fsk.ny,fsk.nx]),
                 vmin=mp_min,vmax=mp_max,
                 origin='lower', interpolation='nearest')
    im.cmap.set_under('#FFFFFF')
    cbar=plt.colorbar(im,orientation='horizontal')
    cbar.ax.set_xlabel(label,fontsize=14)
    yticks = np.array([-6., -5., -4.])
    ax.coords[1].set_ticks(yticks * u.deg)
    ax.coords[1].set_major_formatter('dd')
    ax.set_xlabel('R.A.', fontsize=14)
    ax.set_ylabel('Dec.', fontsize=14)
    if savename is not None:
        plt.savefig(savename,bbox_inches='tight')
plot_syst('Airmass','airmass_maps.fits','g',savename='../doc/Paper/figures/syst_airmass_g.pdf')
plot_syst('Airmass','airmass_maps.fits','i',savename='../doc/Paper/figures/syst_airmass_i.pdf')
plot_syst('Airmass','airmass_maps.fits','y',savename='../doc/Paper/figures/syst_airmass_y.pdf')
plot_syst('Seeing','seeing_maps.fits','g',units='arcsec.',savename='../doc/Paper/figures/syst_seeing_g.pdf',scale=0.17)
plot_syst('Seeing','seeing_maps.fits','i',units='arcsec.',savename='../doc/Paper/figures/syst_seeing_i.pdf',scale=0.17)
plot_syst('Seeing','seeing_maps.fits','y',units='arcsec.',savename='../doc/Paper/figures/syst_seeing_y.pdf',scale=0.17)
plot_syst('Sky count level','skylevel_maps.fits','g',savename='../doc/Paper/figures/syst_skylevel_g.pdf')
plot_syst('Sky count level','skylevel_maps.fits','i',savename='../doc/Paper/figures/syst_skylevel_i.pdf')
plot_syst('Sky count level','skylevel_maps.fits','y',savename='../doc/Paper/figures/syst_skylevel_y.pdf')
plot_syst('Depth','depth_map.fits',None,units='mag.',savename='../doc/Paper/figures/syst_depth.pdf')
plot_syst('Dust','dust_map.fits','g',units='mag.',savename='../doc/Paper/figures/syst_dust.pdf')
plot_syst('Star count','star_map.fits',None,savename='../doc/Paper/figures/syst_star.pdf')

plt.show()
