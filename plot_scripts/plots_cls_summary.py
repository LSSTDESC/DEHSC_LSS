import numpy as np
import matplotlib.pyplot as plt
import sacc
from scipy.stats import chi2 as chi2f

prefix="/global/cscratch1/sd/damonge/HSC_ceci/"

def read_cls(field,
             msk_name='sirius',
             directory='CovAna_NoiAna_MskSirius_ClFit_Dpj1_DpjBands1'):
    """Reads all power spectra for a given field"""
    file_prefix=prefix+"WIDE_"+field+"_"+msk_name+'_i24p5_out/'+directory+'/'
    
    s=sacc.SACC.loadFromHDF(file_prefix+"power_spectra_wdpj.sacc")
    sn=sacc.SACC.loadFromHDF(file_prefix+"noi_bias.sacc")
    sn.precision=s.precision
    sb=sacc.SACC.loadFromHDF(file_prefix+"dpj_bias.sacc")
    sb.precision=s.precision

    return {'data':s,'noise':sn,'dpj_bias':sb}

# Read power spectra for each field and coadd
field_names=['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']
s_fields={}
for f in field_names:
    print(f)
    s_fields[f]=read_cls(f)
s_all={typ:sacc.coadd([s_fields[f][typ] for f in field_names],mode='area')
       for typ in ['data','noise','dpj_bias']}

lmaxs=np.array([2170,2515,3185,4017])

fig,axes=plt.subplots(4,4,figsize=(15,10),sharex=True,sharey=True)
plt.subplots_adjust(hspace=0,wspace=0)

nbins=4
for b1 in range(nbins) :
    for b2 in range(nbins) :
        if b2<b1 :
            axes[b1,b2].axis('off')

for t1,t2,typ,ells,ndx in s_all['data'].sortTracers():
    lmax=min(lmaxs[t1],lmaxs[t2])
    ax=axes[t1,t2]
    cl=s_all['data'].mean.vector[ndx]
    nl=s_all['noise'].mean.vector[ndx]
    dl=s_all['dpj_bias'].mean.vector[ndx]
    el=np.sqrt(np.fabs(np.diag(s_all['data'].precision.getCovarianceMatrix()[ndx,:][:,ndx])))
    if t1==0 and t2==0:
        label_c='Signal'
        label_n='Noise'
        label_d='Deprojection bias'
    else:
        label_c=None
        label_n=None
        label_d=None
    l_factor=1.0
    ax.errorbar(ells,l_factor*(cl-nl),yerr=l_factor*el,fmt='r.',label=label_c)
    eb=ax.errorbar(ells,-l_factor*(cl-nl),yerr=l_factor*el,fmt='r.')
    eb[-1][0].set_linestyle('--')
    ax.plot(ells,l_factor*nl,'-',c='#AAAAAA',label=label_n)
    ax.plot(ells,l_factor*dl,'-',lw=1,c='#0000CC',label=label_d)
    ax.plot(ells,-l_factor*dl,'--',lw=1,c='#0000CC')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.text(0.82,0.9,'$(%d,%d)$'%(t1+1,t2+1),transform=ax.transAxes,fontsize=12)
    ax.set_xlim([110,18000])
    ax.axvspan(lmax,ax.get_xlim()[1],color='grey',alpha=0.2)
    if t1==t2 :
        ax.set_xlabel("$\\ell$",fontsize=16)
        ax.set_ylabel("$C_\\ell$",fontsize=16)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
    else :
        ax.tick_params(axis='both', which='both',length=0)

fig.legend(loc=(0.15,0.15),frameon=False,ncol=1,fontsize=16)
plt.savefig('../doc/Paper/figures/cls_summary.pdf',
            bbox_inches='tight')
plt.show()
