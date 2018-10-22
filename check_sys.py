from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy.stats import binned_statistic
import flatmaps as fm
import os
from astropy.io import fits
from optparse import OptionParser

prefix_data = '/global/cscratch1/sd/damonge/HSC/HSC_processed'
# Define options
parser = OptionParser()

parser.add_option('--input-prefix', dest='prefix_in', default=prefix_data, type=str,
                 help='Input prefix')
parser.add_option('--output-prefix', dest='prefix_out', default='', type=str,
                 help='Output prefix')
parser.add_option('--nsys-bins', dest='nbins_sys', default=7, type=int,
                 help='Number of bins for the systematic analyses')
parser.add_option('--map-path', dest='fname_maps', type=str,
                 help='Path to maps to analyze')
parser.add_option('--depth-cut', dest='depth_cut', default=24.5, type=float,
                 help='Depth cut')
parser.add_option('--mask-threshold', dest='mask_thr', default=0.5, type=float,
                 help='Minimum area fraction of a given pixel to be considered in the analysis')
o, args = parser.parse_args()

def stats_on_sysmap(sys_map, mask, data_map, nbins, bintype='percentiles', perc0=30, reweight=True):
    """ Auxiliary routine that reads a galaxy density map and an observing condition
    density map and computes mean and standard deviation of the galaxy overdensity map
    as a function of the observing condition overdensity map 
    
    Args:
    -----
    
    sys_map: (flatmap) density map of the observing condition that we want to analyze.
    mask: (flatmap) Survey/region mask.
    data_map: (flatmap) Galaxy density map.
    nbins: (int) Number of bins in which to analyze sys_map
    bintype: (str) Binning type: `percentiles` uses the percentiles of sys_map, `equal` makes equal spaced
    bins in sys_map and `log` makes logarithmic spaced bins (default `percentiles`).
    perc0: (float) In case of using the percentiles as binning scheme, starting percentile value (default=30).
    reweight: (bool) If True, reweights the galaxy density map by the mask value (default=True).

    Returns:
    --------
   
    bin_centers: (float) Centers of the bins where we analyze sys_map.
    mean: (float) Mean value of data_map at bin_centers.
    err: (float) Uncertainty on the mean value of data_map at bin_centers.
    """
    if bintype not in ['percentiles', 'equal', 'log']:
        raise ValueError('Only `percentiles`, `equal` or `log` bintypes allowed.')
    mean = np.zeros(nbins)
    err = np.zeros(nbins)
    bin_centers = np.zeros(nbins)
    binary_mask = mask > 0
    if reweight: 
        data_map[binary_mask] = data_map[binary_mask]/mask[binary_mask]
    if bintype=='percentiles':
        for i in range(nbins):
            bin_centers[i] = np.percentile(sys_map[binary_mask], perc0+(100.-perc0)/nbins*i)
            sys_mask = sys_map[binary_mask] < bin_centers[i]
            mean[i] = np.mean(data_map[binary_mask][sys_mask])
            err[i] = np.std(data_map[binary_mask][sys_mask])/np.sqrt(np.count_nonzero(sys_mask))
    if bintype in ['equal','log']:
        if bintype=='equal':
            nbins=nbins
        else:
            min_sys = np.nanmin(sys_map[binary_mask])
            max_sys = np.nanmax(sys_map[binary_mask])
            if min_sys < 0:
                nbins = np.logspace(np.log10(min_sys-min_sys), np.log10(max_sys-min_sys), nbins+1)
            else:
                nbins = np.logspace(np.log10(min_sys), np.log10(max_sys), nbins+1)
        mean, bin_edges, _  = binned_statistic(sys_map[mask], data_map[binary_mask], statistic='mean', bins=nbins)
        std, bin_edges, _  = binned_statistic(sys_map[mask], data_map[binary_mask], statistic='std', bins=nbins)
        counts, bin_edges, _  = binned_statistic(sys_map[mask], data_map[binary_mask], statistic='count', bins=nbins)
        bin_centers = 0.5*bin_edges[1:]+0.5*bin_edges[:-1]
        err = std/np.sqrt(counts)
    return bin_centers, mean, err
     
def check_sys(data_hdu, path_sys, mask, nbins, **kwargs):
    """ Routine to check the evolution of the mean
    density as a function of different potential
    sources of systematic biases/uncertanty on 
    galaxy clustering
  
    Args:
    -----
    
    data_hdu: (HDU) HDU containing the data that we want to analyze.

    path_sys: (str) Path to the flatmap of the contaminant(s).

    mask: (flatmap) Mask of the region/survey to use.

    **kwargs: (dict) arguments to pass to `stats_on_sysmap`

    Outputs:
    --------

    xsys: Values of the potential source of systematics in each bin
    mean: Mean galaxy density in each bin
    err: Error on the mean density in each bin
    """

    fmi, s_map = fm.read_flat_map(path_sys, i_map=-1)
    fmd, data_map = fm.read_flat_map(None, hdu=data_hdu) 
    mean = []
    err = []
    bin_centers = []
    for sys_map in s_map:
        aux_centers, aux_mean, aux_err = stats_on_sysmap(sys_map, mask, data_map, nbins, **kwargs) 
        mean.append(aux_mean)
        bin_centers.append(aux_centers)
        err.append(aux_err)
    return bin_centers, mean, err

# Set up

band = ['u','g','r','i','z']
cont_maps = ['oc_airmass','oc_ccdtemp','oc_ellipt','oc_exptime','oc_nvisit', \
    'oc_seeing', 'oc_sigma_sky', 'oc_skylevel','syst_dust','syst_nstar_i24.50']
xlabels= ['Airmass', r'CCD Temperature [$^{\circ}$C]', 'PSF Ellipticity', \
    'Exposure Time [s]', 'Number of visits', 'Seeing [pixels]', 'Sky noise [ADU]', \
    'Sky level [ADU]', 'E(B-V)', 'Stars per pixel'] 
data_hdus = fits.open(o.fname_maps)
if len(data_hdus)%2!=0:
    raise ValueError("Input file should have two HDUs per map")
nbins = len(data_hdus)//2

print("Reading mask")
#Create depth-based mask
fsk,mp_depth=fm.read_flat_map(o.prefix_in+"_10s_depth_mean_fluxerr.fits",2)
mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=o.depth_cut]=1

#Read masked fraction
fskb,mskfrac=fm.read_flat_map(o.prefix_in+'_MaskedFraction.fits',i_map=0)
fm.compare_infos(fsk,fskb)

#Create BO-based mask
msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>o.mask_thr]=1
msk_t=msk_bo*msk_depth*mskfrac
for ibin in range(nbins):
    data_hdu = data_hdus[2*ibin]
    for j, cm in enumerate(cont_maps):
        path_sys = o.prefix_in+'_%s.fits' %(cm)
        xsys, mean_sys, std_sys = check_sys(data_hdu, path_sys, msk_t, o.nbins_sys) 
        plt.figure()
        if len(xsys)>1:
            for i in range(len(xsys)):
                plt.errorbar(xsys[i], mean_sys[i], std_sys[i], fmt='o', label='%s-band' %band[i], fillstyle='none')
        else:
           plt.errorbar(xsys[0], mean_sys[0], std_sys[0], fmt='o', label=xlabels[j], fillstyle='none')
        plt.ylabel(r'$\bar{n}$ [galaxies/pixel]', fontsize=16)
        plt.xlabel(xlabels[j], fontsize=16)
        plt.grid()
        plt.legend(loc='best')
        mean_sys_all = np.concatenate(mean_sys).ravel()
        ymin = 0.9*np.percentile(mean_sys_all[~np.isnan(mean_sys_all)], 5)
        ymax = 1.1*np.percentile(mean_sys_all[~np.isnan(mean_sys_all)], 95)
        plt.ylim(ymin, ymax)
        if len(xsys)>1:
            sname = cm.split('_')[-1]+'_bin_%d.pdf' % ibin
        else:
            sname = 'nstars_bin_%d.pdf' % ibin
        plt.savefig(os.path.join(o.prefix_out, sname))
        plt.close() 
