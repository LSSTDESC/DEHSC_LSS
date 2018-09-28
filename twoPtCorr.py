from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table
from createMaps import createCountsMap
from optparse import OptionParser
import flatmaps as fm
import sacc
import sys
import time
import os
import treecorr

# We define some convenience functions below

def runTreeCorr(data_ra, data_dec, random_ra, random_dec, # all in degrees
                minSep, maxSep, nBins, data2_ra=None, data2_dec=None,
                random2_ra=None, random2_dec=None):

    startTime= time.time()
    dataCatalog= treecorr.Catalog(ra= data_ra, dec=  data_dec,
                                  ra_units='degrees', dec_units='degrees')
    randCatalog = treecorr.Catalog(ra= random_ra, dec= random_dec,
                                   ra_units='degrees', dec_units='degrees')
    if data2_ra is not None:
        data2Catalog = treecorr.Catalog(ra = data2_ra, dec = data2_dec,
                                        ra_units='degrees', dec_units='degrees')
        rand2Catalog = treecorr.Catalog(ra = rnd2_ra, dec = rnd2_dec,
                                        ra_units='degrees', dec_units='degrees')

    D1D1 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
    D1R1 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
    R1R1 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')

    # use LS
    D1D1.process(dataCatalog, randCatalog)
    R1R1.process(randCatalog, randCatalog)
    D1R1.process(dataCatalog, randCatalog)
    w1, vx1 = D1D1.calculateXi(R1R1, D1R1)
    theta = np.exp(DD.logr)

    if data2_ra is not None:
        D1D2 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        R1R2 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        R2D2 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        D1D2 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        D2D2 = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        R2R2 = treecorr.NNCorrleation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units= 'degrees')
        D1D2.process(dataCatalog, data2Catalog)
        R1R2.process(randCatalog, rand2Catalog)
        R2D2.process(data2Catalog, rand2Catalog)
        R2R2.process(rand2Catalog, rand2Catalog)
        D2D2.process(data2Catalog, data2Catalog)
        w12, vx12 = D1D2.calculateXi(R1R2, D1R1, R2D2)
        w2, vx2 = D2D2.calculateXi(R2R2, D2R2)
    else:
        w12 = None
        w2 = None
 
    print 'Time taken: %s (s)'%(time.time()-startTime)

    return theta, w1, w12, w2

def create_randoms(mask, nrnd):
    """ Function to generate randoms"""
def compute_cov(ra, dec, ra_rnd, dec_rnd, minSep, maxSep, nBins, ra2=None, dec2=None, ra2_rnd=None, dec2_rnd=None, njk=100):
    """ Function to compute the JK covariance"""
    import kmeans_radec
    from kmeans_radec import KMeans
    X = np.zeros((len(ra),2))
    X[:,0] = ra
    X[:,1] = dec
    X2 = np.zeros((len(ra_rnd),2))
    X2[:,0] = ra_rnd
    X2[:,1] = dec_rnd
    cen_guess = np.zeros((njk,2))
    xmin = np.nanmin(ra)
    xmax = np.nanmax(ra)
    ymin = np.nanmin(dec)
    ymax = np.nanmin(dec)
    # Initial guess for the centers
    cen_guess[:,0] = np.linspace(xmin,xmax,njk)
    cen_guess[:,1] = np.linspace(ymin,ymax,njk)
    km = KMeans(cen_guess,tol=1e-3)
    km.run(X, maxiter=100)
    labels = km.find_nearest(X)
    labels_rnd = km.find_nearest(X2)
    def run_corr(label):
        theta, w, dw = runTreeCorr(,dec,ra_rnd,o.max_sep,o.nbins)
        return theta, w, dw
    print('Starting covariance calculation')
    t0 = time()
    pool = multiprocessing.Pool(processes=o.n_proc)
    # We are going to perform a take one out JK
    param = np.unique(labels)
    th_list, w_list, dw_list = pool.map(run_corr,param)
    w_list = np.array(w_list)
    cov = np.covariance(w_list,rowvar=False)


##################################################

# Setting options and path


prefix_data='/global/cscratch1/sd/damonge/HSC/'
def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
#Options

parser.add_option('--input-prefix', dest='prefix_in', default='NONE', type=str,
                  help='Input prefix. '+
                  'Systematics maps will be input-prefix + _<syst_name>.fits. '+
                  'Mask will be input-prefix + _MaskedFraction.fits')
parser.add_option('--input-catalog', dest='fname_catalog', default='NONE', type=str,
                  help='Path to input catalog file')
parser.add_option('--input-random', dest='fname_randoms', default='NONE', type=str,
                  help='Path to input randoms file')
parser.add_option('--theta-min', dest='theta_min', default=1e-3, type=double,
                  help='Minimum pair separation to consider for the correlation function [degrees]')
parser.add_option('--theta-max', dest='theta_max', default=5, type=double,
                  help='Maximum pair separation to consider for the correlation function [degrees]')
parser.add_option('--nbins-theta', dest='nbins_theta', default=15, type=int,
                  help='Number of angular bins in which to compute the correlation function')
parser.add_option('--output-file', dest='fname_out',default=None,type=str,
                  help='Output file name. Output will be in SACC format.')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
                  help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--depth-cut', dest='depth_cut', default=24.5, type=float,
                  help='Minimum depth to consider in your footprint')
parser.add_option('--masking-threshold',dest='mask_thr',default=0.5, type=float,
                  help='Will discard all pixel with a masked fraction larger than this.')
parser.add_option('--hsc-field',dest='hsc_field',default='HSC_WIDE',type=str,
                  help="HSC field used here (just for labelling purposes)")
parser.add_option('--cont-depth',dest='cont_depth',default=False,action='store_true',
                  help='Marginalize over depth map template')
parser.add_option('--cont-dust',dest='cont_dust',default=False,action='store_true',
                  help='Marginalize over dust template')
parser.add_option('--cont-dust-bands',dest='cont_dust_bands',default=False,
                  action='store_true',help='Marginalize over dust template in all bands')
parser.add_option('--cont-stars',dest='cont_stars',default=False,action='store_true',
                  help='Marginalize over stars template')
parser.add_option('--cont-oc',dest='cont_oc',type='string',
                  action='callback',callback=opt_callback,
                  help='If you want to marginalize over observing condition maps, '+
                  'list here all observing conditions you want to include')
parser.add_option('--cont-oc-bands',dest='cont_oc_bands',default=False,action='store_true',
                  help='Marginalize over observing contition templates in all bands.')
parser.add_option('--pz-type',dest='pz_type',default='nnpz',type=str,
                  help='Photo-z to use')
parser.add_option('--pz-mark',dest='pz_mark',default='best',type=str,
                  help='Photo-z summary statistic to use when binning objects')
parser.add_option('--pz-bins',dest='fname_bins',default=None,type=str,
                  help='File containing the redshift bins (format: 1 row per bin, 2 columns: z_ini z_end)')
parser.add_option('--do-cross', dest='do_cross', default=False, action='store_true',
                  help='Compute angular cross-correlations')
parser.add_option('--n-jk', dest='njk', default=100, type=int,
                  help='Number of JK regions to compute the covariances')


######################################################


# Read options
(o, args) = parser.parse_args()

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
msk_t=msk_bo*msk_depth

#Area
area_patch=np.sum(msk_t*mskfrac)*np.radians(fsk.dx)*np.radians(fsk.dy)

#Read contaminants
print("Reading contaminant templates")
temps=[]
def read_map_bands(fname,read_bands) :
  if read_bands :
    i_map=-1
  else :
    i_map=['g','r','i','z','y'].index(o.band)
  fskb,temp=fm.read_flat_map(fname,i_map=i_map)
  fm.compare_infos(fsk,fskb)
  if i_map!=-1 :
    temp=[temp]
    
  return temp
# 1- Depth
if o.cont_depth :
  temps.append(mp_depth)
# 2- Dust
if o.cont_dust :
  temp=read_map_bands(o.prefix_in+'_syst_dust.fits',o.cont_dust_bands)
  for t in temp :
    temps.append(t)
# 3- Stars
if o.cont_stars :
  fskb,temp=fm.read_flat_map(o.prefix_in+'_syst_nstar_'+o.band+'%.2lf.fits'%o.depth_cut,
                             i_map=0)
  fm.compare_infos(fsk,fskb)
  temps.append(temp)
# 4- Observing conditions
if o.cont_oc is not None :
  oc_all=['airmass','ccdtemp','ellipt','exptime','nvisit','seeing','sigma_sky','skylevel']
  #Include only OCs we know about
  oc_list=[]
  for c in o.cont_oc :
      if c in oc_all :
        oc_list.append(c)
  for c in oc_list :
    temp=read_map_bands(o.prefix_in+'_oc_'+c+'.fits',o.cont_oc_bands)
    for t in temp :
      temps.append(t)
if temps==[] :
  temps=None
else :
  print(" - Will marginalize over a total of %d contaminant templates"%(len(temps)))

# Photo-z stuff

if (o.fname_bins is None) or (not os.path.isfile(o.fname_bins)) :
  raise KeyError("Can't find bins file")

if o.pz_type=='ephor_ab' :
  pz_code='eab'
elif o.pz_type=='frankenz' :
  pz_code='frz'
elif o.pz_type=='nnpz' :
  pz_code='nnz'
else :
  raise KeyError("Photo-z method "+o.pz_type+" unavailable. Choose ephor_ab, frankenz or nnpz")

if o.pz_mark  not in ['best','mean','mode','mc'] :
  raise KeyError("Photo-z mark "+o.pz_mark+" unavailable. Choose between best, mean, mode and mc")

column_mark='pz_'+o.pz_mark+'_'+pz_code
column_pdfs='pz_mc_'+pz_code

column_mark_rnd = column_mark # We'll assume that the randoms also have some sort of photo-z associated?

zi_arr,zf_arr=np.loadtxt(o.fname_bins,unpack=True,ndmin=2)
nbins=len(zi_arr)

# Read catalog produced by process.py

cat = Table.read(o.fname_catalog)

# Read randoms or create them

if o.fname_randoms=='NONE':
    ra_rnd, dec_rnd, z_rnd = create_randoms(mask)
else:
    rnd =  Table.read(o.fname_randoms)
    ra_rnd = rnd['ra']
    dec_rnd = rnd['dec']
    z_rnd = rnd[column_mark_rnd]

# Compute auto and x-corrs on the data points
w1_all, w12_all, w2_all = []
if o.do_cross:
    for i in range(nbins):
        m1 = (cat[column_mark]>zi_arr[i]) & (cat[column_mark]<=zf_arr[i])
        m1r = (z_rnd>zi_arr[i]) & (z_rnd<=zf_arr[i])
        ra1 = cat['ra'][m1]
        dec1 = cat['dec'][m1]
        for j in range(i,nbins):
            m2 = (cat[column_mark]>zi_arr[j]) & (cat[column_mark]<=zf_arr[j])
            m2r = (z_rnd>zi_arr[j]) & (z_rnd<zf_arr[j])
            ra2 = cat['ra'][m2]
            dec2 = cat['dec'][m2]
            theta, w1_aux, w12_aux, w2_aux = runTreeCorr(ra1, dec1, ra_rnd[m1r], dec_rnd[m1r], o.theta_min, o.theta_max, o.nbins_theta, 
                                            ra2, dec2, ra_rnd[m2r], dec_rnd[m2r])
            w1_all.append(w1_aux)
            w2_all.append(w2_aux)
            w12_all.append(w12_aux)
else:                 
    for bin_low, bin_high in zip(zi_arr,zf_arr):
        ra1 = cat['ra'][(cat[column_mark]>bin_low) & (cat[column_mark]<=bin_high)]
        dec1 = cat['dec'][(cato[column_mark]>bin_low) & (cat[column_mark]<=bin_high)]
        m1r = (z_rnd>bin_low) & (z_rnd<=bin_high)
        theta, w1_aux, w12_aux, w2_aux = runTreeCorr(ra1, dec1, ra_rnd[m1r], dec_rnd[m1r], o.theta_min, o.theta_max, o.nbins_theta)
        w1_all.append(w1_aux)
        w2_all.append(w2_aux)
        w12_all.append(w12_aux)

