import numpy as np
import matplotlib.pyplot as plt
import sacc
from scipy.stats import chi2 as chi2f
import formatting

prefix="/global/cscratch1/sd/damonge/HSC_ceci/"

def read_cls(field,
             wdpj=True,
             msk_name='sirius',
             directory='CovAna_NoiAna_MskSirius_ClFit_Dpj0_DpjBands1'):
    file_prefix=prefix+"WIDE_"+field+"_"+msk_name+'_i24p5_out/'+directory+'/'
    if msk_name=='arcturus':
        file_prefix='../data_replotting/arcturus/'+field+'/'+directory+'/'
    else:
        file_prefix='../data_replotting/'+field+'/'+directory+'/'

    if wdpj:
        s=sacc.SACC.loadFromHDF(file_prefix+"power_spectra_wdpj.sacc")
    else:
        s=sacc.SACC.loadFromHDF(file_prefix+"power_spectra_wodpj.sacc")
    sn=sacc.SACC.loadFromHDF(file_prefix+"noi_bias.sacc")
    sn.precision=s.precision
    sb=sacc.SACC.loadFromHDF(file_prefix+"dpj_bias.sacc")
    sb.precision=s.precision

    return {'data':s,'noise':sn,'dpj_bias':sb}

def read_cls_coadd(msk_name,directory,wdpj=True):
    field_names=['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']
    s_fields={}
    print(msk_name,directory)
    for f in field_names:
        print("   %s "%f)
        s_fields[f]=read_cls(f,msk_name=msk_name,directory=directory,wdpj=wdpj)
    s_all={typ:sacc.coadd([s_fields[f][typ] for f in field_names],mode='area')
           for typ in ['data','noise','dpj_bias']}
    return s_all

s_fid=read_cls_coadd('sirius','CovAna_NoiAna_MskSirius_ClFit_Dpj0_DpjBands1_newDepth',wdpj=True)
s_ndp=read_cls_coadd('sirius','CovAna_NoiAna_MskSirius_ClFit_Dpj0_DpjBands1_newDepth',wdpj=False)
s_sys=read_cls_coadd('sirius','CovAna_NoiAna_MskSiriusSyst_ClFit_Dpj0_DpjBands1_newDepth',wdpj=True)
s_msk=read_cls_coadd('arcturus','CovAna_NoiAna_MskArcturus_ClFit_Dpj0_DpjBands1_newDepth',wdpj=True)

lmaxs=np.array([2170,2515,3185,4017])

fig,axes=plt.subplots(4,4,figsize=(15,10),sharex=True,sharey=True)
plt.subplots_adjust(hspace=0,wspace=0)

nbins=4
for b1 in range(nbins) :
    for b2 in range(nbins) :
        if b2<b1 :
            axes[b2,b1].axis('off')

for t1,t2,typ,ells,ndx in s_fid['data'].sortTracers():
    lmax=min(lmaxs[t1],lmaxs[t2])
    mask=ells<lmax
    ls=ells[mask]
    ndxs=ndx[mask]
    ax=axes[t2,t1]
    c_fid=s_fid['data'].mean.vector[ndxs]
    c_ndp=s_ndp['data'].mean.vector[ndxs]
    c_sys=s_sys['data'].mean.vector[ndxs]
    c_msk=s_msk['data'].mean.vector[ndxs]
    e_fid=np.sqrt(np.fabs(np.diag(s_fid['data'].precision.getCovarianceMatrix()[ndxs,:][:,ndxs])))
    e_ndp=np.sqrt(np.fabs(np.diag(s_ndp['data'].precision.getCovarianceMatrix()[ndxs,:][:,ndxs])))
    e_sys=np.sqrt(np.fabs(np.diag(s_sys['data'].precision.getCovarianceMatrix()[ndxs,:][:,ndxs])))
    e_msk=np.sqrt(np.fabs(np.diag(s_msk['data'].precision.getCovarianceMatrix()[ndxs,:][:,ndxs])))
    print("%d-%d"%(t1+1,t2+1))
    ax.plot([2,20000],[0,0],'k--')
    if t1==0 and t2==0:
        label_ndp='No deprojection'
        label_sys='Contaminant mask'
        label_msk='Arcturus mask'
    else:
        label_ndp=None
        label_sys=None
        label_msk=None
    ax.errorbar(ls*(0.95+0.1*(0+0.5)/3.),(c_ndp-c_fid)/e_ndp,
                yerr=e_ndp/e_ndp,
                fmt='.',c=formatting.cmap1[0],label=label_ndp,
                markersize=5, elinewidth=1.5, capthick=1.5, capsize=3.5)
                
    ax.errorbar(ls*(0.95+0.1*(1+0.5)/3.),(c_sys-c_fid)/e_sys,
                yerr=e_sys/e_sys,
                fmt='.',c=formatting.cmap1[1],label=label_sys,
                markersize=5, elinewidth=1.5, capthick=1.5, capsize=3.5)
    ax.errorbar(ls*(0.95+0.1*(2+0.5)/3.),(c_msk-c_fid)/e_sys,
                yerr=e_sys/e_sys,
                fmt='.',c=formatting.cmap1[2],label=label_msk,
                markersize=5, elinewidth=1.5, capthick=1.5, capsize=3.5)
    ax.set_xscale('log')
    ax.set_ylim([-3.9,3.9])
    ax.set_xlim([90,4000])
    ax.text(0.05,0.9,'$(%d,%d)$'%(t1+1,t2+1),transform=ax.transAxes,fontsize=12)
    #if t1==t2 :
    ax.set_xlabel("$\\ell$",fontsize=18)
    ax.set_ylabel("$\\Delta C_\\ell/\\sigma(C_\\ell)$",fontsize=18)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    ##else :
    #    ax.tick_params(axis='both', which='both',length=0)

fig.legend(loc=(0.75,0.85),frameon=False,ncol=1,fontsize=20)
plt.savefig('../doc/Paper/figures/cls_systematics.pdf',
            bbox_inches='tight')

plt.show()
