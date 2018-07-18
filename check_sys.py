from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy.stats import binned_statistic
import flatmaps as fm
import os
import astropy.table

def check_sys(path_data, path_sys, mask_path, nbins=20):
    """ Routine to check the evolution of the mean
    density as a function of different potential
    sources of systematic biases/uncertanty on 
    galaxy clustering
  
    Args:
    -----
    
    path_data: Path to the catalog containing the clean data,
    it should contain RA and DEC.

    path_sys: Path to the flatmap of the contaminant(s)

    mask_path: Path to file containing the mask

    nbins: Number of bins in which to perform the analysis
    (20 by default)

    Outputs:
    --------

    xsys: Values of the potential source of systematics in each bin
    mean: Mean galaxy density in each bin
    err: Error on the mean density in each bin
    """

    fmi, s_map = fm.read_flat_map(path_sys, i_map=-1)
    fmm, mask = fm.read_flat_map(mask_path)
    binary_mask = mask>0
    data = astropy.table.Table.read(path_data)
    data_map = np.bincount(fmi.pos2pix(data['ra'],data['dec']),minlength=fmi.npix)
    mean = []
    err = []
    bin_centers = []
    for sys_map in s_map:
        aux_mean, bin_edges, _  = binned_statistic(sys_map[binary_mask], data_map[binary_mask]/mask[binary_mask], statistic='mean', bins=nbins)
        aux_std, bin_edges, _  = binned_statistic(sys_map[binary_mask], data_map[binary_mask]/mask[binary_mask], statistic='std', bins=nbins)
        aux_counts, bin_edges, _  = binned_statistic(sys_map[binary_mask], data_map[binary_mask]/mask[binary_mask], statistic='count', bins=nbins)
        mean.append(aux_mean)
        bin_centers.append(0.5*(bin_edges[1:]+bin_edges[:-1]))
        err.append(1.0*aux_std/np.sqrt(aux_counts))
    return bin_centers, mean, err

path_base = '/global/cscratch1/sd/damonge/HSC/HSC_processed'
field = 'WIDE_GAMA15H'
path_data = os.path.join(path_base,field,'%s_Catalog_i24.50.fits' %field)
path_mask = os.path.join(path_base,field,'%s_MaskedFraction.fits' %field)
band = ['u','g','r','i','z']
cont_maps = ['oc_airmass','oc_ccdtemp','oc_ellipt','oc_exptime','oc_nvisit', \
    'oc_seeing', 'oc_sigma_sky', 'oc_skylevel','syst_dust','syst_nstar_i24.50']
xlabels= ['Airmass', r'CCD Temperature [$^{\circ}$C]', 'PSF Ellipticity', \
    'Exposure Time [s]', 'Number of visits', 'Seeing [pixels]', 'Sky noise [ADU]', \
    'Sky level [ADU]', 'E(B-V)', 'Stars per pixel'] 
for j, cm in enumerate(cont_maps):
    path_sys = os.path.join(path_base,field,'%s_%s.fits' %(field,cm))
    xsys, mean_sys, std_sys = check_sys(path_data, path_sys, path_mask)
    plt.figure()
    if len(xsys)>1:
        for i in range(len(xsys)):
            plt.errorbar(xsys[i],mean_sys[i],std_sys[i],fmt='o',label='%s-band' %band[i],fillstyle='none')
    else:
       plt.errorbar(xsys[0],mean_sys[0],std_sys[0],fmt='o',label=xlabels[j],fillstyle='none')
    plt.ylabel(r'$\bar{n}$ [galaxies/pixel]',fontsize=16)
    plt.xlabel(xlabels[j],fontsize=16)
    plt.grid()
    plt.legend(loc='best')
    plt.ylim(5,20)
    if len(xsys)>1:
        sname = cm.split('_')[-1]+'.pdf'
    else:
        sname = 'nstars.pdf'
    plt.savefig(sname) 
