from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.table
from createMaps import createMeanStdMaps, createCountsMap
import estDepth
from optparse import OptionParser
import flatmaps as fm
import sys
import intermediates
from twoPtCorr import *
import time
def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))
parser = OptionParser()
# Options
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
parser.add_option('--harmonic', dest='harmonic', default=False, action='store_true',
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
parser.add_option('--use-fitsio',dest='use_fitsio',default=False,action='store_true',
    help='Use fitsio to open fits files')
parser.add_option('--n-rnd',dest='n_rnd',default=10000000,type=int,
    help='Number of random particles to compute 2pt correlations')
parser.add_option('--min-sep',dest='min_sep',default=0.005,type=float,
    help='Minimum separation angle to compute the correlation function (in deg)')
parser.add_option('--max-sep',dest='max_sep',default=10,type=float,
    help='Maximum separation angle to compute the correlation function (in deg)')
parser.add_option('--num-bins',dest='num_bins',default=25,type=int,
    help='Number of bins to compute the correlation function')
parser.add_option('--do-jk',dest='do_jk',default=False,action='store_true',
    help='Perform jack-knife')
parser.add_option('--num-reg',dest='num_reg',default=100,type=int,
    help='Number of regions for JK')
parser.add_option('--num-proc',dest='num_proc',default=4,type=int,
    help='Number of parallel JK processes')
parser.add_option('--output-cov',dest='fname_cov',default=None,
    help='Name of the output covariance')
parser.add_option("--templates", dest="templates_fname", default="none",type="string",
    help="Templates to subtract from power-spectrum",action="callback", callback=opt_callback)

####
# Read options
(o, args) = parser.parse_args()
# Read catalog (with fitsio or astropy)
if o.use_fitsio==True:
    import fitsio
    try:
        cat = fitsio.read(o.fname_in)
    except:
        raise TypeError('%s is does not contain a FITS table/specify input file' % o.fname_in)
        sys.exit(1)
if o.use_fitsio==False:
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
# Define bounds for the skymap, we add two pixels on the edges 
minx = np.nanmin(cat['ra'])-2*o.res
miny = np.nanmin(cat['dec'])+2*o.res
maxx = np.nanmax(cat['ra'])-2*o.res
maxy = np.nanmax(cat['dec'])+2*o.res
print('Map edges: ', minx, maxx, miny, maxy)
# Generate flat sky map that will contain our galaxies
flatSkyGrid = fm.FlatMapInfo([minx,maxx],[miny,maxy],dx=o.res,dy=o.res)
# Clean the catalog to estimate the depths
sel = np.ones(len(cat),dtype=bool)
print('Read ', len(cat), ' objects')
if o.use_fitsio==False:
    for tab_key in cat.keys():
        sel = sel & (np.isnan(cat[tab_key])==False)
else:
    for tab_key in cat.dtype.fields.keys():
        sel = sel & (np.isnan(cat[tab_key])==False)
print('Selected ', np.count_nonzero(sel), ' objects')
# Check if we want to save the depth maps
if o.dirname_out_depth is not None:
   save_depth=True
else:
   save_depth=False
# Compute SNR
snr = cat['%scmodel_flux'%o.band][sel]/cat['%scmodel_flux_err'%o.band][sel]
# Create depth maps
if o.method==0:
    depth, depth_std = estDepth.dr1paper_method(cat['ra'][sel],cat['dec'][sel],o.band,mags=cat['%scmodel_mag'%o.band][sel], \
    snr=snr,flatSkyGrid=flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)
if o.method==1:
    depth, depth_std = estDepth.desc_method(cat['ra'][sel],cat['dec'][sel],o.band,mags=cat['%scmodel_mag'%o.band][sel], \
    snr=snr,flatSkyGrid=flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)
if o.method==2:
    depth, depth_std = estDepth.flux_err(cat['ra'][sel],cat['dec'][sel],cat['%scmodel_flux_err'%o.band][sel], \
    o.band,flatSkyGrid,SNRthreshold=o.min_snr,plotMaps=o.gen_plot, \
    saveMaps=save_depth,otuputDir=o.dirname_out_depth)
if o.method==3:
    depth, depth_std = estDepth.depth_map_meanSNRrange(cat['ra'][sel],cat['dec'][sel],o.band,mags=cat['%scmodel_mag'%o.band][sel], \
    snr=snr,flatSkyGrid=flatSkyGrid,SNRthreshold=o.min_snr, \
    plotMaps=o.gen_plot,saveMaps=save_depth,outputDir=o.dirname_out_depth)
# Now we have the depth maps we can select the footprint
fp_mask = depth >= o.depth_cut
if o.gen_plot:
    flatSkyGrid.view_map(fp_mask,posColorbar= True,
                          title= 'Mask %s-band %sigma depth'%(o.band, o.min_snr),
                          xlabel='ra', ylabel='dec')
    plt.show()
# Make the galaxy selection
good_gals = (cat['%scmodel_mag'%o.band] >= o.min_mag) & (cat['%scmodel_mag'%o.band] <= o.max_mag) & (cat['iclassification_extendedness'] == 1)
pixnums = flatSkyGrid.pos2pix(cat['ra'][good_gals],cat['dec'][good_gals])
in_footprint = np.in1d(pixnums,np.where(fp_mask)[0])
ra_data = cat['ra'][good_gals][in_footprint]
dec_data  = cat['dec'][good_gals][in_footprint]
# Make random catalog
if o.harmonic==False:
    randCatalog = intermediates.getRandomCatalog(flatSkyGrid,fp_mask,minRA=minx, \
        maxRA=maxx, minDec=miny, maxDec=maxy, nData = o.n_rnd, plotMap=o.gen_plot)
    theta, w, dw = runTreeCorr(data_ra = ra_data, data_dec = dec_data, \
        random_ra = randCatalog['ra'], random_dec = randCatalog['dec'], \
        minSep = o.min_sep, maxSep = o.max_sep, nBins= o.num_bins)
    if o.gen_plot:
        plt.plot(theta,w,'o')
        plt.xlabel(r'$\theta$ [deg]')
        plt.ylabel(r'$w(\theta)$')
        plt.loglog()
        plt.show()
    tab_2pt = astropy.table.Table([theta,w,dw],names=('theta','w','dw'))
    tab_2pt.write(o.fname_out_2pt,overwrite=True)
    if o.do_jk:
        import multiprocessing
        labels, labels_rnd = get_labels(ra_data,dec_data,randCatalog['ra'],randCatalog['dec'],o.num_reg)
        
        def run_corr(label):
            theta, w, dw = runTreeCorr(data_ra=ra_data[labels!=label],data_dec=dec_data[labels!=label], \
                random_ra = randCatalog['ra'][labels_rnd!=label], random_dec = randCatalog['dec'][labels_rnd!=label], minSep=o.min_sep, \
              maxSep=o.max_sep, nBins=o.num_bins)
            return w
        print('Starting covariance calculation')
        t0 = time.time()
        pool = multiprocessing.Pool(processes=o.num_proc)
        # We are going to perform a take one out JK
        param = np.unique(labels)
        print('Going to perform ',len(param), ' JK regions')
        w_list = pool.map(run_corr,param)
        w_list = np.array(w_list)
        cov = np.cov(w_list,rowvar=False)
        if o.gen_plot:
            plt.figure()
            plt.imshow(cov)
            plt.show()
        t1 = time.time()
        print('Done with covariances. Ellapsed time %f (s)' % (t1-t0))
        np.save(o.fname_cov,cov)
    # Calculate the cross-correlations with systematics
    #if "none" not in o.templates_fname:
    #    temps=[]
    #    for tname in o.templates_fname:
    #        temps.append(createSysMap(flatSkyGrid,hp.read_map(tname),fp_mask,plotMaps=o.gen_plot))
if o.harmonic:
    mp = createCountsMap(ra_data, dec_data, flatSkyGrid, returnMap= True)
    dmap = np.zeros(len(mp))
    dmap[mp!=0]=mp[mp!=0]/np.mean(mp[mp!=0])-1.
    cl, lbpw, wsp = flatSkyGrid.compute_power_spectrum(dmap*fp_mask,fp_mask)
    ells = np.mean(lbpw,axis=0)
    tab_2pt = astropy.table.Table([ells,cl],names=('l','cl'))
    tab_2pt.write(o.fname_out_2pt,overwrite=True)
    if o.do_jk:
        import multiprocessing
        labels = get_labels_one(ra_data,dec_data,o.num_reg)
        def run_cl(label):
            mp = createCountsMap(ra_data[labels!=label], dec_data[labels!=label], flatSkyGrid, returnMap= True)
            dmap = np.zeros(len(mp))
            dmap[mp!=0]=mp[mp!=0]/np.mean(mp[mp!=0])-1.
            cl, lbpw, wsp = flatSkyGrid.compute_power_spectrum(dmap*fp_mask,fp_mask)
            return cl
        print('Starting covariance calculation')
        t0 = time.time()
        pool = multiprocessing.Pool(processes=o.num_proc)
        # We are going to perform a take one out JK
        param = np.unique(labels)
        print('Going to perform ',len(param), ' JK regions')
        cl_list = pool.map(run_cl,param)
        cl_list = np.array(cl_list)
        cov = np.cov(cl_list,rowvar=False)
        if o.gen_plot:
            plt.figure()
            plt.imshow(cov)
            plt.show()
        t1 = time.time()
        print('Done with covariances. Ellapsed time %f (s)' % (t1-t0))
        np.save(o.fname_cov,cov) 
