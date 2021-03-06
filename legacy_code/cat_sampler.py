from __future__ import print_function
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from createMaps import createCountsMap
from optparse import OptionParser
import flatmaps as fm
import sys
import time
import os

def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
#Options
parser.add_option('--input-prefix', dest='prefix_in', default='NONE', type=str,
                  help='Input prefix. The input catalog will be searched for as input-prefix + _Catalog_<band><limit>.fits.')
parser.add_option('--output-file', dest='fname_out',default=None,type=str,
                  help='Output file name. If None, I\'ll use input-prefix + _bins_ + fname_bins + .fits')
parser.add_option('--no-bo-cut',dest='no_bo_cut',default=False,action='store_true',
                  help='Do not remove objects within bright-object mask')
parser.add_option('--pz-type',dest='pz_type',default='nnpz',type=str,
                  help='Photo-z to use')
parser.add_option('--pz-mark',dest='pz_mark',default='best',type=str,
                  help='Photo-z summary statistic to use when binning objects')
parser.add_option('--use-pdf',dest='use_pdf',default=False,
                  help='Whether to stack photo-z pdfs to generate N(z)')
parser.add_option('--use-cosmos',dest='use_cosmos',default=False,
                  help='Whether to use the COSMOS reweighting to generate N(z)')
parser.add_option('--pz-bins',dest='fname_bins',default=None,type=str,
                  help='File containing the redshift bins (format: 1 row per bin, 2 columns: z_ini z_end)')
parser.add_option('--nz-bins',dest='nz_bin_num',default=200,type=int,
                  help='Number of bins to use in the output N(z)\'s')
parser.add_option('--nz-max',dest='nz_bin_max',default=4.0,type=float,
                  help='Maximum redshift to use for output N(z)\'s')
parser.add_option('--map-sample',dest='map_sample',default=None,type=str,
                  help='Sample map used to determine the pixelization that will be used. If None I\'ll try to find the masked fraction map')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
                  help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--depth-cut', dest='depth_cut', default=24.5, type=float,
                  help='Minimum depth to consider in your footprint')
parser.add_option('--cosmos-weights-file',dest='cosmos_weights_file',default='/global/cscratch1/sd/damonge/HSC/HSC_processed/COSMOS_HSC_WEIGHTS.fits',type=str,
                  help='Photo-z summary statistic to use when binning objects')
parser.add_option('--bo-mask-type', dest='mask_type', default='sirius', type=str,
                  help='Bright object mask (arcturus or sirius)')

####
# Read options
(o, args) = parser.parse_args()

prefix_in_use=o.prefix_in
if o.mask_type=='arcturus' :
  prefix_in_use+='_marct'

if o.use_cosmos and o.use_pdf :
  raise KeyError("Can't do both pdf stacking and COSMOS 30-band")
  exit(1)

fname_cat=prefix_in_use+'_Catalog_'+o.band+'%.2lf.fits'%o.depth_cut
if not os.path.isfile(fname_cat) :
  raise KeyError("File "+fname_cat+" doesn't exist")

if o.map_sample is None :
  o.map_sample=prefix_in_use+'_MaskedFraction.fits'
if not os.path.isfile(o.map_sample) :
  raise KeyError("File "+o.map_sample+" doesn't exist")

if (o.fname_bins is None) or (not os.path.isfile(o.fname_bins)) :
  raise KeyError("Can't find bins file")

if o.fname_out is None :
  o.fname_out=prefix_in_use+'_bins_'+o.fname_bins+'.fits'

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
column_zmc='pz_mc_'+pz_code

#Read catalog
cat=fits.open(fname_cat)[1].data
if not o.no_bo_cut :
  if o.mask_type=='arcturus' :
    msk=cat['mask_Arcturus'].astype(bool)
  else :
    msk=np.logical_not(cat['iflags_pixel_bright_object_center'])
    msk*=np.logical_not(cat['iflags_pixel_bright_object_any'])
  cat=cat[msk]

if o.use_pdf:
  # Read in pdfs and bins
  pdf_file = o.prefix_in+'_pdfs_'+o.pz_type+'.fits'
  pdfs = fits.open(pdf_file)[1].data['pdf']
  bins = fits.open(pdf_file)[2].data['bins']
  pdfs = pdfs[msk]

if o.use_cosmos:
  # Read in weights and bins
  weights_file = fits.open(o.cosmos_weights_file)[1].data

#Read map information
fsk,mpdum=fm.read_flat_map(o.map_sample,0)

#Read bins
zi_arr,zf_arr=np.loadtxt(o.fname_bins,unpack=True,ndmin=2)
nbins=len(zi_arr)

#Iterate through bins
maps=[]
nzs=[]
for zi,zf in zip(zi_arr,zf_arr) :
  msk=(cat[column_mark]<=zf) & (cat[column_mark]>zi)
  subcat=cat[msk]
  if o.use_pdf:
    binpdfs = pdfs[msk] # The pdfs which have a redshift in the correct bin
    bz = np.linspace(0,o.nz_bin_max,o.nz_bin_num+1)
    hz = []
    for x in xrange(len(bz) - 1):
      # interpolate at the edges of the bins and then integrate
      redshift_subset = [max(bz[x], bins[0])] + list(bins[(bins > bz[x]) & (bins < bz[x+1])]) + [bz[x+1]] # The interpolation x-axis
      interp_pdfs = interp1d(bins, binpdfs)(redshift_subset) # Get the PDF values at the interpolation points
      pdf_area = np.trapz(interp_pdfs, x = redshift_subset, axis = 1)
      hz.append(np.nansum(pdf_area))
    hz = np.array(hz)
  elif o.use_cosmos:
    msk_cosmos=(weights_file['hsc_'+column_mark]<=zf) & (weights_file['hsc_'+column_mark]>zi)
    binweights=weights_file[msk_cosmos]['weight']
    bincosmosz=weights_file[msk_cosmos]['cosmos_photoz']
    hz,bz=np.histogram(bincosmosz,bins=o.nz_bin_num,range=[0.,o.nz_bin_max],
                       weights=binweights)
    hnz,bnz=np.histogram(bincosmosz,bins=o.nz_bin_num,range=[0.,o.nz_bin_max])
    ehz=np.zeros(len(hnz)); ehz[hnz>0]=(hz[hnz>0]+0.)/np.sqrt(hnz[hnz>0]+0.)
  else:
    zmcs=subcat[column_zmc]
    hz,bz=np.histogram(zmcs,bins=o.nz_bin_num,range=[0.,o.nz_bin_max])
  nmap=createCountsMap(subcat['ra'],subcat['dec'],fsk)
  if o.use_cosmos :
    nzs.append([bz[:-1],bz[1:],hz+0.,ehz])
  else :
    nzs.append([bz[:-1],bz[1:],hz+0.])
  maps.append(nmap)
nzs=np.array(nzs)
maps=np.array(maps)

#Save maps and N(z)s
if len(maps[0])!=fsk.npix :
  raise ValueError("Map doesn't conform to this pixelization")

header=fsk.wcs.to_header()
hdus=[]
for im,m in enumerate(maps) :
  #Map
  head=header.copy()
  head['DESCR']=('Ngal, bin %d'%(im+1),'Description')
  if im==0 :
    hdu=fits.PrimaryHDU(data=m.reshape([fsk.ny,fsk.nx]),header=head)
  else :
    hdu=fits.ImageHDU(data=m.reshape([fsk.ny,fsk.nx]),header=head)
  hdus.append(hdu)

  #Nz
  cols=[fits.Column(name='z_i',array=nzs[im,0,:],format='E'),
        fits.Column(name='z_f',array=nzs[im,1,:],format='E'),
        fits.Column(name='n_z',array=nzs[im,2,:],format='E')]
  if o.use_cosmos :
    cols.append(fits.Column(name='en_z',array=nzs[im,3,:],format='E'))
    
  hdus.append(fits.BinTableHDU.from_columns(cols))
hdulist=fits.HDUList(hdus)
hdulist.writeto(o.fname_out,overwrite=True)
