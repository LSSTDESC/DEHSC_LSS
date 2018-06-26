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
import os

def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
# Options
parser.add_option('--input-file', dest='fname_in', default='NONE', type=str,
                  help='Path to input FITS table')
parser.add_option('--output-file',dest='fname_out', default=None, type=str,
                  help='Path to output directory')

####
# Read options
(o, args) = parser.parse_args()
# Read catalog (with fitsio or astropy)
print("Reading")
try:
  data = Table.read(o.fname_in)
except:
  raise TypeError('%s is does not contain a FITS table/specify input file' % o.fname_in)

####
# Clean nulls and nans
print("Basic cleanup")
sel=np.ones(len(data),dtype=bool)
names=[n for n in data.keys()]
isnull_names=[]
for key in data.keys() :
  if key.__contains__('isnull') :
    sel[data[key]]=0
    isnull_names.append(key)
data.remove_columns(isnull_names)
data.remove_rows(~sel) #np.where(~sel)[0])#[sel]

# Remove useless columns
fields_save=['frame_id','frame_num','exp_id','exp_num','ccd_id','ccd','ccdname','pointing',
             'ccdtemp','ra','decl','equinox','ra2000','decl2000',
             'llcra','llcdecl','ulcra','ulcdecl','urcra','urcdecl','lrcra','lrcdecl',
             'pa','insrot','mjd','azimuth','elevation','airmass','exptime',
             'gain1','gain2','gain3','gain4',
             'skylevel','sigma_sky','seeing','ellipt','ell_pa','wcs_nobj','wcs_rms',
             'fluxmag0','fluxmag0err','zeropt','zeropt_err','nobj_bright','flag_auto','filter']
to_remove=[]
for key in data.keys() :
  if key not in fields_save :
    to_remove.append(key)
data.remove_columns(to_remove)

print(len(data),len(data.keys()))

# Write out
hdr=fits.Header()
hdr['PREPROC']='yes'
prm_hdu=fits.PrimaryHDU(header=hdr)
dat_hdu=fits.table_to_hdu(data)
hdul=fits.HDUList([prm_hdu,dat_hdu])
hdul.writeto(o.fname_out,clobber=True)
