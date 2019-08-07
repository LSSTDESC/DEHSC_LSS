import flatmaps as fm
import numpy as np
import matplotlib.pyplot as plt

for field in ['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']:
    fsk,msk=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+field+"_sirius_i24p5_out/masked_fraction.fits")
    fig=plt.figure()
    ax=fig.add_subplot(111,projection=fsk.wcs)
    ax.set_title(field,fontsize=14)
    im=ax.imshow(msk.reshape([fsk.ny,fsk.nx]),vmin=0,vmax=1,
                 origin='lower', interpolation='nearest')
    ax.set_xlabel('R.A.', fontsize=14)
    ax.set_ylabel('Dec.', fontsize=14)
    plt.savefig("../doc/Paper/figures/mask_"+field+".pdf",bbox_inches='tight')
plt.show()
