from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.table
from createMaps import createMeanStdMaps
import estDepth
from optparse import OptionParser
import flatmaps as fm
import sys
parser = OptionParser()

parser.add_option('--input-catalog', dest='fname_in', default=None, type=str,
    help='Path to input catalog (FITS table)')
parser.add_option('--min-mag', dest='min_mag', default=15, type=float,
    help='Minimum magnitude to consider in the analysis')
parser.add_option('--max-mag', dest='max_mag', default=25.0, type=float,
    help='Maximum magnitude to consider in the analysis')
parser.add_option('--min-snr', dest='min_snr', default=10.0, type=float,
    help='SNR threshold used to compute the depth maps')
parser.add_option('--output-mask', dest='fname_out_mask', default=None, type=str,
    help='Path to output mask')
parser.add_option('--output-depth-dir',dest='dirname_out_depth', default=None, type=str,
    help='Path to output depth map') 
parser.add_option('--output-catalog', dest='fname_out_cat', default=None, type=str,
    help='Path to output (clean) catalog')
parser.add_option('--fourier', dest='fourier', default=False, action='store_true',
    help='If selected compute the power-spectrum')
parser.add_option('--output-2pt', dest='fname_out_2pt', default=None, type=str,
    help='Path to output 2pt results')
parser.add_option('--plot', dest='gen_plot', default=False, action='store_true',
    help='If selected show and save some plots')
parser.add_option('--depth-cut', dest='depth_cut', default=25.0, type=float,
    help='Minimum depth to consider in your footprint')
parser.add_option('--resolution', dest='res', default=0.0285, type=float,
    help='Map/mask resolution (in degrees)')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
    help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--method', dest='method', default='0', type=int,
    help='Method to construct depth maps: 0-> DR1-like, 1-DESC, 2-Mean SNR, 3-Mean flux')
# Read options
(o, args) = parser.parse_args()
# Read catalog
try:
    cat = astropy.table.Table.read(o.fname_in)
except:
    raise TypeError('%s is does not contain a FITS table/specify input file' % o.fname_in)
    sys.exit(1)
bands = ['g','r','i','z','y']
# Check if the band is available in the catalog
if o.band not in bands:
   print('Selected band not available, select g, r, i, z or y')
   sys.exit(1)
# Define bounds for the skymap
minx = np.min(cat['ra'])
miny = np.min(cat['dec'])
maxx = np.max(cat['ra'])
maxy = np.max(cat['dec'])
# Generate flat sky map that will contain our galaxies
flatSkyGrid = fm.FlatMapInfo([minx,maxx],[miny,maxy],dx=o.res,dy=o.res)
# Clean the catalog to estimate the depths
sel = np.ones(len(cat),dtype=bool)
print('Read ', len(cat), ' objects')
for tab_key in cat.keys():
    sel = sel & (np.isnan(cat[tab_key])==False)
print('Selected ', np.count_nonzero(sel), ' objects')
# Check if we want to save the depth maps
if o.dirname_out_depth is not None:
   save_depth=True
else:
   save_depth=False
# Create depth maps
if o.method==0:
    estDepth.dr1paper_method(cat['ra'][sel],cat['dec'][sel],o.band,cat['%scmodel_mag'%o.band][sel], \
    cat['%scmodel-SNR'%o.band][sel],flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)
if o.method==1:
    estDepth.desc_method(cat['ra'][sel],cat['dec'][sel],o.band,cat['%scmodel_mag'%o.band][sel], \
    cat['%scmodel-SNR'%o.band][sel],flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)
if o.method==2:
    estDepth.flux_err(cat['ra'][sel],cat['dec'][sel],cat['%scmodel_flux_err'%o.band][sel], \
    o.band,flatSkyGrid,SNRthreshold=o.min_snr,plotMaps=o.gen_plot, \
    saveMaps=save_depth,otuputDir=o.dirname_out_depth)
if o.method==3:
    estDepth.depth_map_meanSNRrange(cat['ra'][sel],cat['dec'][sel],o.band,cat['%scmodel_mag'%o.band][sel], \
    cat['%scmodel-SNR'%o.band][sel],flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)

