from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from createMaps import createCountsMap
from optparse import OptionParser
import flatmaps as fm
import sys
import time
import os

prefix_data='/global/cscratch1/sd/damonge/HSC/'
def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
#Options
parser.add_option('--inout-prefix', dest='prefix_io', default='NONE', type=str,
                  help='Input prefix. The input catalog will be searched for as inout-prefix + _Catalog_<band><limit>.fits.')
parser.add_option('--no-bo-cut',dest='no_bo_cut',default=False,action='store_true',
                  help='Remove objects within bright-object mask')
parser.add_option('--pz-type',dest='pz_type',default='nnpz',type=str,
                  help='Photo-z to use')
parser.add_option('--pz-mark',dest='pz_mark',default='best',type=str,
                  help='Photo-z summary statistic to use when binning objects')
parser.add_option('--pz-bins',dest='fname_bins',default=None,type=str,
                  help='File containing the redshift bins (format: 1 row per bin, 2 columns: z_ini z_end)')
parser.add_option('--map-sample',dest='map_sample',default=None,type=str,
                  help='Sample map used to determine the pixelization that will be used. If None I\'ll try to find the masked fraction map')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
                  help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--depth-cut', dest='depth_cut', default=24.5, type=float,
                  help='Minimum depth to consider in your footprint')

####
# Read options
(o, args) = parser.parse_args()

fname_field=o.prefix_io+'_Catalog_'+o.band+'%.2lf.fits'%o.depth_cut
if not os.path.isfile(fmame_field) :
    raise KeyError("File "+fname_field+" doesn't exist")

if o.map_sample is None :
    o.map_sample=o.prefix_io+'_MaskedFraction.fits'
if os.path.isfile(o.map_sampe) :
    raise KeyError("File "+o.map_sample+" doesn't exist")

if (o.fname_bins is None) or (not os.path.isfile(o.fname_bins)) :
    raise KeyError("Can't fine bins file")

if o.pz_type=='ephor_ab' :
    pz_code='eab'
elif o.pz_type='frankenz' :
    pz_code='frz'
elif o.pz_type='nnpz' :
    pz_code='nnz'
else :
    raise KeyError("Photo-z method "+o.pz_type+" unavailable. Choose ephor_ab, frankenz or nnpz")

if o.pz_mark  not in ['best','mean','mode','mc'] :
    raise KeyError("Photo-z mark "+o.pz_mark+" unavailable. Choose between best, mean, mode and mc")

column_mark='pz_'+o.pz_mark+'_'+pz_code

print(column_mark)

print(o.prefix_io)
print(o.no_bo_cut)
print(o.pz_type)
print(o.pz_mark)
print(o.fname_bins)
print(o.map_sample)
