from __future__ import print_function
import numpy as np
import os
import flatmaps as fm
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from astropy.io import fits
import intermediates
from twoPtCorr import runTreeCorr
import pickle
from plots import plot_wthetas
import os

HSCprocessedDir="/global/cscratch1/sd/damonge/HSC/HSC_processed/"
SNRthreshold= 10
magCut= 24.5
band= 'i'
fields=['WIDE_GAMA09H','WIDE_GAMA15H', 'WIDE_HECTOMAP', 'WIDE_VVDS',
        'WIDE_WIDE12H','WIDE_AEGIS','WIDE_XMMLSS']
#fields=['WIDE_GAMA15H']

fstar= 0.01  # reading off left panel in Fig 17 in https://arxiv.org/pdf/1705.06766.pdf
minSep= 0.001; maxSep= 5; nBins= 25

outDir= '/global/cscratch1/sd/awan/HSCOutput/'

tbins, w, w_sig= {}, {}, {}

for field in fields:
    print('\n***\nLooking at %s'%field)
    # HSCprocessedDir has files in folder with name= field. access fits, npzs only.
    HSCdatapath= '%s/%s/'%(HSCprocessedDir, field)
    HSCFiles= [f for f in os.listdir(HSCdatapath) if (f.__contains__('fits') or f.__contains__('npz'))]
    HSCFiles= [HSCdatapath+f for f in HSCFiles]

    # read in the relevant datafiles
    # read in depth map
    for filename in  [f for f in HSCFiles if f.__contains__('%ss_depth_mean'%SNRthreshold)]:
        print('Reading ', filename)
        flatSkyGrid, depthMap= fm.read_flat_map(filename)
    # read in dust map
    #for filename in  [f for f in HSCFiles if f.__contains__('syst_dust')]:
    #    print('reading ', filename)
    #    newGrid, dustMap= fm.read_flat_map(filename)
    #    if newGrid.get_dims()!=flatSkyGrid.get_dims() :
    #        raise ValueError("Inconsistent fskys")
    # read in star map
    for filename in  [f for f in HSCFiles if f.__contains__('_syst_nstar_i%.2lf.npz'%magCut)]:
        print('Reading ', filename)
        newGrid, starMap= fm.read_flat_map(filename)
        if newGrid.get_dims()!=flatSkyGrid.get_dims() :
            raise ValueError("Inconsistent fskys")
            
    # clean up depth map as in study_power_spectra.py
    depthMap[np.isnan(depthMap)]=0
    depthMap[depthMap>40]=0
    
    # -------
    # construct galaxy catalog
    # read in the catalog
    filename= [f for f in HSCFiles if f.__contains__('_Catalog_i%.2lf.fits'%magCut) ][0]
    print('Reading ', filename)
    catalog= fits.open(filename)[1].data 
    print('\nNobj in the read-in galaxy catalog: ', np.shape(catalog))
    
    # apply thr bright objects maps; following study_power_spectra.py
    mask_BO= ~(catalog['iflags_pixel_bright_object_center']*catalog['iflags_pixel_bright_object_any'])
    catalog= catalog[mask_BO]
    print('Nobj after BO mask: ', np.shape(catalog))

    # find the mag-limited galaxy sample
    ind= np.where((catalog['%scmodel_mag'%band]<magCut) & (catalog['iclassification_extendedness']==1))[0]
    galSample_ra, galSample_dec= catalog['ra_r'][ind], catalog['dec_r'][ind]
    print('%s galaxies to in the mag-limited sample.'%len(galSample_ra))
    
    # -------
    # create a star catalog from the mag limited star map
    fg_u, starMap_ug= flatSkyGrid.u_grade(starMap, x_fac= 5)  # increase the pixelization by 5x
    pixsWithStars= np.where(starMap_ug>0)[0]
    print('\nConsidering pixels with max %s stars'%max(starMap_ug[pixsWithStars]))
    starSample_ra, starSample_dec= fg_u.pix2pos(pixsWithStars)
    print('%s stars to consider.\n'%len(starSample_ra))
    
    # -------
    # create random depth map
    mask= intermediates.getMask(depthMap, band, magCut, plotMap= False,
                                flatSkyGrid= flatSkyGrid, title= '')

    # create random catalog
    nData= max(len(galSample_ra), len(starSample_ra))   # max output star, gal sample. probably too much.
    randCatalog= intermediates.getRandomCatalog(flatSkyGrid, mask, 
                                                minRA= min(catalog['ra_r']), maxRA= max(catalog['ra_r']),
                                                minDec=  min(catalog['dec_r']), maxDec=  max(catalog['dec_r']),
                                                nData= nData, plotMap= False)
    if (len(randCatalog)>0):
        # -----
        # run tree corr. nn correlation.
        print('\nRunning TreeCorr ... ')
        theta, wtheta, wtheta_sig= {}, {}, {}

        tag= 'gal_noCorr'
        print('\nRunning for', tag)
        theta[tag], wtheta[tag], wtheta_sig[tag]= runTreeCorr(data_ra= galSample_ra, data_dec= galSample_dec,
                                                              random_ra= randCatalog['ra'], random_dec= randCatalog['dec'],
                                                              minSep= minSep, maxSep= maxSep, nBins= nBins)
        tag= 'stars'
        print('\nRunning for', tag)
        theta[tag], wtheta[tag], wtheta_sig[tag]= runTreeCorr(data_ra= starSample_ra, data_dec= starSample_dec,
                                                              random_ra= randCatalog['ra'], random_dec= randCatalog['dec'],
                                                              minSep= minSep, maxSep= maxSep, nBins= nBins)
        tag= 'gal_starCorrected'
        print('\nRunning for', tag)
        # correction based on eq. 26 in https://arxiv.org/abs/1507.05360
        theta[tag]= theta['gal_noCorr']
        wtheta[tag]= (1+fstar)**2*(wtheta['gal_noCorr']-fstar**2*wtheta['stars'])-fstar**4
        wtheta_sig[tag]= np.sqrt((1+fstar)**2*(wtheta_sig['gal_noCorr']**2-fstar**2*wtheta_sig['stars']**2))


        # -----
        # save the corrs
        filename= 'correlations_%s.pickle'%field

        current= os.getcwd()
        os.chdir(outDir)
        with open(filename, 'wb') as handle:
            pickle.dump({'theta': theta, 'wtheta': wtheta,
                        'wtheta_sig': wtheta_sig},
                        handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.chdir(current)
        print('Saved ', filename)
        # -----
        # plot
        plot_wthetas(theta, wtheta, wtheta_sig, title= '%s: %s<%s'%(field, band, magCut),
                     savePlot= True, outDir= outDir, filename= 'correlations_%s_interms.png'%field)

        # save the final-correction one for comparison across fields
        tag= 'gal_starCorrected'
        tbins[field], w[field], w_sig[field]= theta[tag], wtheta[tag], wtheta_sig[tag]
    else:
        print('*** !!!! SOMETHING IS WRONG WITH %s.'%field)
        
# plot all final ones
plot_wthetas(tbins, w, w_sig, title= '%s<%s'%(band, magCut),
            savePlot= True, outDir= outDir, filename= 'correlations_allcorrected.png')
print('All done.')
