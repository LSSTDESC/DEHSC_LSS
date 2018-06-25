from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table,vstack
from astropy import wcs
from createMaps import createMeanStdMaps, createCountsMap, removeDisconnected
from flatMask import createMask
import estDepth
from optparse import OptionParser
import flatmaps as fm
import sys
import time
import rotate as rtt
import os

bands = ['g','r','i','z','y']
depth_methods = ['fluxerr'] #dr1','desc','fluxerr']
prefix_data='/global/cscratch1/sd/damonge/HSC/'
def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
# Options
parser.add_option('--input-field', dest='field_in', default='NONE', type=str,
                  help='Path to input catalog (FITS table)')
parser.add_option('--min-snr', dest='min_snr', default=10.0, type=float,
                  help='SNR threshold used to compute the depth maps')
parser.add_option('--output-prefix',dest='out_prefix', default=None, type=str,
                  help='Path to output directory')
parser.add_option('--save-depth-maps', dest='sv_depth', default=False, action='store_true',
                  help='If selected save depth maps')
parser.add_option('--save-masks', dest='sv_mask', default=False, action='store_true',
                  help='If selected save masks')
parser.add_option('--save-systematics', dest='sv_syst', default=False, action='store_true',
                  help='If selected save systematics maps')
parser.add_option('--gen-plots', dest='gen_plot', default=False, action='store_true',
                  help='If selected create and save some plots')
parser.add_option('--show-plots', dest='show_plot', default=False, action='store_true',
                  help='If selected show some plots')
parser.add_option('--depth-cut', dest='depth_cut', default=25.0, type=float,
                  help='Minimum depth to consider in your footprint')
parser.add_option('--resolution', dest='res', default=0.0285, type=float,
                  help='Map/mask resolution (in degrees)')
parser.add_option('--resolution-bo', dest='res_bo', default=0.003, type=float,
                  help='B.O. mask resolution (in degrees)')
parser.add_option('--field-padding', dest='pad', default=-1, type=float,
                  help='Edge padding (in degrees)')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
                  help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--depth-method', dest='depth_method', default='0', type=int,
                  help='Method to construct depth maps: 0-> DR1-like, 1-DESC, 2-Flux-error')

####
# Read options
(o, args) = parser.parse_args()
# Read catalog (with fitsio or astropy)
print("Reading")
nparts=0
while os.path.isfile(prefix_data+'HSC_'+o.field_in+'_part%d_forced.fits'%(nparts+1)) :
  nparts+=1
if nparts==0 :
  #Single-part file
  fname_in=prefix_data+'HSC_'+o.field_in+'_forced.fits'
  try:
    cat = Table.read(fname_in)
  except:
    raise TypeError('%s is does not contain a FITS table/specify input file' % fname_in)
    sys.exit(1)
else :
  #Parted file
  for ipart in np.arange(nparts)+1 :
    print(" Reading part %d"%ipart)
    fname_in=prefix_data+'HSC_'+o.field_in+'_part%d_forced.fits'%ipart
    try:
      c=Table.read(fname_in)
    except:
      raise TypeError('%s is does not contain a FITS table/specify input file' % fname_in)
      sys.exit(1)
    print(len(c))
    if ipart==1 :
      cat=c
    else :
      cat=vstack([cat,c],join_type='exact')
print(len(cat))

# Check if the band is available in the catalog
if o.band not in bands:
   print('Selected band not available, select g, r, i, z or y')
   sys.exit(1)

####
# Clean nulls and nans
print("Basic cleanup")
sel=np.ones(len(cat),dtype=bool)
names=[n for n in cat.keys()]
isnull_names=[]
for key in cat.keys() :
  if key.__contains__('isnull') :
    sel[cat[key]]=0
    isnull_names.append(key)
  else :
    sel[np.isnan(cat[key])]=0
print("Will drop %d rows"%(len(sel)-np.sum(sel)))
cat.remove_columns(isnull_names)
cat.remove_rows(~sel) #np.where(~sel)[0])#[sel]

####
# Generate flat-sky information object
fsk=fm.FlatMapInfo.from_coords(cat['ra'],cat['dec'],o.res,pad=o.pad/o.res)

####
# Generate systematics maps
syst={}
systdesc={}
# 1- Dust
dustmaps=[]
dustdesc=[]
for b in bands :
  m,s=createMeanStdMaps(cat['ra'],cat['dec'],cat['a_'+b],fsk,nan_outside=False)
  if (b==o.band) and o.gen_plot :
    fsk.view_map(m,posColorbar= True,title= 'A_%s'%b,
                 xlabel='ra', ylabel='dec',colorMin=np.amin(m[m>0]),
                 fnameOut=o.out_prefix+'_a_%s.png'%b)
  dustmaps.append(m)
  dustdesc.append('Dust, '+b+'-band')
syst['dust']=np.array(dustmaps)
systdesc['dust']=np.array(dustdesc)

# 2- Nstar
#    This needs to be done for stars passing the same cuts as the sample (except for the s/g separator)
sel_maglim=np.ones(len(cat),dtype=bool); sel_maglim[cat['%scmodel_mag'%o.band]-cat['a_%s'%o.band]>o.depth_cut]=0
sel_stars=np.ones(len(cat),dtype=bool);  sel_stars[cat['iclassification_extendedness']>0.99]=0
sel_gals =np.ones(len(cat),dtype=bool);  sel_gals[cat['iclassification_extendedness']<0.99]=0
mstar=createCountsMap(cat['ra'][sel_maglim*sel_stars],cat['dec'][sel_maglim*sel_stars],fsk)+0.
if o.gen_plot :
  fsk.view_map(mstar,posColorbar=True,title='N_star',xlabel='ra',ylabel='dec',
               fnameOut=o.out_prefix+'_nstar_'+o.band+'%.2lf.png'%(o.depth_cut))
syst['nstar_'+o.band+'%.2lf'%(o.depth_cut)]=mstar
systdesc['nstar_'+o.band+'%.2lf'%(o.depth_cut)]='Stars, '+o.band+'<%.2lf'%(o.depth_cut)
if o.sv_syst :
   for k in syst.keys() :
     fsk.write_flat_map(o.out_prefix+'_syst_'+k+'.fits',syst[k],descript=systdesc[k])

####
# Generate bright-object mask
#Binary BO mask
mask_bo,fsg=createMask(cat['ra'],cat['dec'],
                       [cat['iflags_pixel_bright_object_center'],
                        cat['iflags_pixel_bright_object_any']],
                       fsk,o.res_bo)
if o.gen_plot :
  fsg.view_map(mask_bo,posColorbar= True,title= 'Bright-object mask',
               xlabel='ra', ylabel='dec',colorMin=0,colorMax=1,
               fnameOut=o.out_prefix+'_BOMask.png')
if o.sv_mask :
  fsg.write_flat_map(o.out_prefix+'_BOMask.fits',mask_bo,descript='Bright-object mask')

#Masked fraction
masked=np.ones(len(cat))
masked*=np.logical_not(cat['iflags_pixel_bright_object_center'])
masked*=np.logical_not(cat['iflags_pixel_bright_object_center'])
masked_fraction,s=createMeanStdMaps(cat['ra'],cat['dec'],masked,fsk)
masked_fraction_cont=removeDisconnected(masked_fraction,fsk)
if o.gen_plot :
  fsk.view_map(masked_fraction_cont,posColorbar=True,title='Masked fraction',
               xlabel='ra',ylabel='dec',colorMin=0,colorMax=1,
               fnameOut=o.out_prefix+'_MaskedFraction.png')
if o.sv_mask :
  fsk.write_flat_map(o.out_prefix+'_MaskedFraction.fits',masked_fraction_cont,descript='Masked fraction')

####
# Compute SNR and depth maps
print("Creating depth maps")
snrs={};
depths={};
for i,b in enumerate(bands) :
  snrs[b]=cat['%scmodel_flux'%b]/cat['%scmodel_flux_err'%b]
  depths[b]={}
  plot_stuff= (b==o.band) and o.gen_plot
  for method in depth_methods :
    print(b,method,plot_stuff)
    depths[b][method]={}
    if method=='fluxerr' :
      depth,depth_std=estDepth.get_depth(method,cat['ra'],cat['dec'],b,arr1=cat['%scmodel_flux_err'%b],
                                         arr2=None,flatSkyGrid=fsk,SNRthreshold=o.min_snr,
                                         plotMaps=plot_stuff,prefixOut=o.out_prefix)
    else :
      depth,depth_std=estDepth.get_depth(method,cat['ra'],cat['dec'],b,arr1=cat['%scmodel_mag'%b],
                                         arr2=snrs[b],flatSkyGrid=fsk,SNRthreshold=o.min_snr,
                                         plotMaps=plot_stuff,prefixOut=o.out_prefix)
    depths[b][method]['depth']=depth; depths[b][method]['depth_std']=depth_std
if o.sv_depth :
  for method in depth_methods :
    mps_m=[]; mps_s=[];
    dcp_m=[]; dcp_s=[];
    for b in bands :
      mps_m.append(depths[b][method]['depth'])
      mps_s.append(depths[b][method]['depth_std'])
      dcp_m.append('%d-s depth, '%o.min_snr+b+' '+method+' mean')
      dcp_s.append('%d-s depth, '%o.min_snr+b+' '+method+' STD')
    fsk.write_flat_map(o.out_prefix+'_%ds_depth_mean_'%(o.min_snr)+method+'.fits',
                       np.array(mps_m),descript=np.array(dcp_m))
    fsk.write_flat_map(o.out_prefix+'_%ds_depth_std_'%(o.min_snr)+method+'.fits',
                       np.array(mps_s),descript=np.array(dcp_s))

####
# Compute depth-based mask
if o.sv_mask :
  print("Depth-based mask")
  dmp=(depths[o.band][depth_methods[o.depth_method]]['depth']).copy();
  #dmp[dmp>100]=-1 #Remove nans
  msk_depth=np.zeros_like(dmp); msk_depth[dmp>=o.depth_cut]=1
  msk_depth=removeDisconnected(msk_depth,fsk)
  name_plot='_DepthMask_'+depth_methods[o.depth_method]+'_'+o.band+'%.2lf'%o.depth_cut
  if o.gen_plot :
    fsk.view_map(msk_depth,posColorbar=True,title='%s>%.2lf mask'%(o.band,o.depth_cut),
                 xlabel='ra',ylabel='dec',colorMin=0,colorMax=1,
                 fnameOut=o.out_prefix+name_plot+'.png')
  fsk.write_flat_map(o.out_prefix+name_plot+'.fits',msk_depth,descript='Depth mask')

####
# Implement final cuts
# - Mag. limit
# - Star-galaxy separator
print("Will lose %d objects to depth and stars"%(np.sum(sel_maglim*sel_gals)))
cat.remove_rows(~(sel_maglim*sel_gals))

####
# Write final catalog
# 1- header
print("Writing output")
hdr=fits.Header()
hdr['BAND']=o.band
hdr['DEPTH']=o.depth_cut
hdr['FIELDN']=o.field_in
prm_hdu=fits.PrimaryHDU(header=hdr)
# 2- Catalog
cat_hdu=fits.table_to_hdu(cat)
# 3- Actual writing
hdul=fits.HDUList([prm_hdu,cat_hdu])
#hdul.writeto(o.out_prefix+'_Catalog_'+o.band+'%.2lf'%(o.depth_cut)+'.fits',overwrite=True)
hdul.writeto(o.out_prefix+'_Catalog_'+o.band+'%.2lf'%(o.depth_cut)+'.fits',clobber=True)

####
# Generate plots and exit
if o.gen_plot and o.show_plot :
  plt.show();
