import numpy as np
import matplotlib.pyplot as plt
import sacc
from scipy.stats import chi2 as chi2f

prefix="/global/cscratch1/sd/damonge/HSC_ceci/"

def read_cls(field,
             msk_name='sirius',
             directory='CovAna_NoiAna_MskSirius_ClFit_Dpj0_DpjBands1'):
    file_prefix=prefix+"WIDE_"+field+"_"+msk_name+'_i24p5_out/'+directory+'/'
    
    s=sacc.SACC.loadFromHDF(file_prefix+"power_spectra_wdpj.sacc")
    sn=sacc.SACC.loadFromHDF(file_prefix+"noi_bias.sacc")
    sn.precision=s.precision
    sb=sacc.SACC.loadFromHDF(file_prefix+"dpj_bias.sacc")
    sb.precision=s.precision

    return {'data':s,'noise':sn,'dpj_bias':sb}

field_names=['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']
colors={'GAMA09H':'#CC0000',
        'GAMA15H':'#CC6600',
        'HECTOMAP':'#CCCC00',
        'VVDS':'#009900',
        'WIDE12H':'#004C99',
        'XMMLSS':'#990099'}

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
    mask=ells<lmax
    ax=axes[t1,t2]
    dmean=s_all['data'].mean.vector[ndx]-s_all['noise'].mean.vector[ndx]
    print("%d-%d"%(t1+1,t2+1))
    ax.plot([2,20000],[0,0],'k--')
    for i_f,f in enumerate(field_names):
        c=colors[f]
        s=s_fields[f]
        d=s['data'].mean.vector[ndx]-s['noise'].mean.vector[ndx]
        cov=s['data'].precision.getCovarianceMatrix()[ndx,:][:,ndx]
        e=np.sqrt(np.fabs(np.diag(cov)))

        res=(d-dmean)[mask]
        icov_red=np.linalg.inv(cov[mask,:][:,mask])
        chi2=np.einsum('i,ij,j',res,icov_red,res)
        ndof=len(res)
        print("  %s : %.3lf %d %.3lf"%(f,chi2,ndof,1-chi2f.cdf(chi2,ndof)))
        if t1==0 and t2==0:
            label=f
        else:
            label=None
        ax.errorbar(ells[mask]*(0.95+0.1*(i_f+0.5)/6.),((d-dmean)/e)[mask],
                    yerr=np.ones_like(e)[mask],
                    fmt='.',c=c,label=label)
    ax.set_xscale('log')
    ax.set_ylim([-3.9,3.9])
    ax.set_xlim([90,4000])
    ax.text(0.05,0.9,'$(%d,%d)$'%(t1+1,t2+1),transform=ax.transAxes,fontsize=12)
    if t1==t2 :
        ax.set_xlabel("$\\ell$",fontsize=16)
        ax.set_ylabel("$\\Delta C_\\ell/\\sigma(C_\\ell)$",fontsize=16)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
    else :
        ax.tick_params(axis='both', which='both',length=0)
        
fig.legend(loc=(0.15,0.15),frameon=False,ncol=2,fontsize=16)
plt.savefig('../doc/Paper/figures/cls_consistency.pdf',
            bbox_inches='tight')
plt.show()
