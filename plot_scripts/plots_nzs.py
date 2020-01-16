import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import formatting

# Bin summary
pz_bins=[0.15,0.50,0.75,1.00,1.50]
ng_arr=[]
nz_arr=[]
nz_arr={'cosmos':[],
        'demp':[],
        'ephor':[],
        'ephor_ab':[],
        'frankenz':[]}
labels={'cosmos':'COSMOS 30-band',
        'demp':'DEmP',
        'ephor':'Ephor',
        'ephor_ab':'Ephor\_AB',
        'frankenz':'FRANKEN-Z'}
pz_codes=sorted(list(nz_arr.keys()))

for fieldname in ['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']:
    predir = "/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+fieldname+"_sirius_i24p5_out/"
    predir = "../data_replotting/"+fieldname+"/"
    mskd=fits.open(predir+"masked_fraction.fits")[0].data
    mskbin=np.zeros_like(mskd); mskbin[mskd>0.5]=1
    fng=fits.open(predir+"ngal_maps.fits")
    zis=fng[1].data['z_i']
    zfs=fng[1].data['z_f']
    zs=0.5*(fng[1].data['z_f']+fng[1].data['z_i'])
    ngs=[]
    for i in range(4):
        ngs.append(np.sum(fng[2*i].data*mskbin))
    ng_arr.append(ngs)
    for c in pz_codes:
        nzs=[]
        for i in range(4):
            nzs.append(fng[2*i+1].data['nz_'+c])
        nz_arr[c].append(nzs)

ng_arr=np.array(ng_arr)
for c in pz_codes:
    nz_arr[c]=np.array(nz_arr[c])

nzs_total={}
for c in pz_codes:
    nzs_total[c]=np.sum(nz_arr[c]*ng_arr[:,:,None],axis=0)

colz=['#CC0000','#0000CC','#FFAB00','#009900']
colz=formatting.cmap1
fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
plt.subplots_adjust(hspace=0,wspace=0)
def plot_in_ax(ax,code):
    ymax=0
    for iz in range(4):
        pz=nzs_total['cosmos'][iz]/np.sum(nzs_total['cosmos'][iz])
        ymax=max(ymax,np.amax(pz))
        ax.step(zfs,pz,c=colz[iz],lw=2)
        pz=nzs_total[code][iz]/np.sum(nzs_total[code][iz])
        ymax=max(ymax,np.amax(pz))
        ax.plot(zs,pz,'-',c=colz[iz],lw=2)
        ax.set_ylim([0,1.05*ymax])
        ax.set_xlim([0,2.4])
        ax.set_yticks([])
        ax.set_xticks([0.0,0.5,1.0,1.5,2.0])
        ax.text(0.95,0.85,labels[code],horizontalalignment='right',
                transform=ax.transAxes,fontsize=12)
plot_in_ax(axes[0][0],'demp')
plot_in_ax(axes[0][1],'ephor')
plot_in_ax(axes[1][0],'ephor_ab')
plot_in_ax(axes[1][1],'frankenz')
axes[0][0].set_ylabel('$p(z)$',fontsize=20)
axes[1][0].set_ylabel('$p(z)$',fontsize=20)
axes[1][0].set_xlabel('$z$',fontsize=20)
axes[1][1].set_xlabel('$z$',fontsize=20)
axes[1][0].tick_params(labelsize="large")
axes[1][1].tick_params(labelsize="large")
plt.savefig('../doc/Paper/figures/nzs.pdf',
            bbox_inches='tight')
plt.show()
