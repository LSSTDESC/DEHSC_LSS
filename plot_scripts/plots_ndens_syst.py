import numpy as np
import matplotlib.pyplot as plt

prefix="/global/cscratch1/sd/damonge/HSC_ceci/WIDE_VVDS_sirius_i24p5_out/systmap_plots/"

# Read representative data
d_airmass=np.load(prefix+"bin3_airmass.npz")
d_dust=np.load(prefix+"bin3_dust.npz")

# Plot
def plot_ndata(x,y,ey,ax,xcut=None):
    # Linear slope
    m = np.sum((x-1)*(y-1)/ey**2)/np.sum((x-1)**2/ey**2)

    # Plot data
    ax.errorbar(x,y,yerr=ey,fmt='ro',ms=4)
    # Plot linear fit
    x0=x[0]-(x[-1]-x[0])*0.05
    xf=x[-1]+(x[-1]-x[0])*0.05
    ax.plot([x0,xf],[m*(x0-1)+1,m*(xf-1)+1],'k-')
    # Plot cutoff
    if xcut:
        y0=np.amin(y-ey)
        yf=np.amax(y+ey)
        
        ax.plot([xcut,xcut],[y0-0.05*(yf-y0),yf+0.05*(yf-y0)],'--',c='#AAAAAA')
        ax.set_ylim([y0-0.05*(yf-y0),yf+0.05*(yf-y0)])

    #for tick in ax.xaxis.get_major_ticks():
    #    tick.label.set_fontsize(14)
    #for tick in ax.yaxis.get_major_ticks():
    #    tick.label.set_fontsize(14)
    ax.tick_params(labelsize="x-large")
    ax.set_ylabel('$R(N_g)$',fontsize=15)
    
fig,axes=plt.subplots(2,1,figsize=(10,10))
# Plot airmass
plot_ndata(d_airmass['x_rescaled'][1],
           d_airmass['mean'][1],
           d_airmass['error'][1],
           axes[0],xcut=1.15)
axes[0].set_xlabel('$R({\\rm Airmass},r)$',fontsize=15)
# Plot dust
plot_ndata(d_dust['x_rescaled'][0],
           d_dust['mean'][0],
           d_dust['error'][0],
           axes[1])
axes[1].set_xlabel('$R({\\rm Dust})$',fontsize=15)
plt.savefig('../doc/Paper/figures/ndens_syst.pdf',
            bbox_inches='tight')
plt.show()
