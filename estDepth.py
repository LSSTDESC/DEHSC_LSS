import numpy as np
from createMaps import createMeanStdMaps

def saveDepthMaps(depth, depth_std, methodName, outputDir, SNRthreshold, band, flatSkyGrid, ):
    if outputDir is None:
        print 'outputDir not specified so maps will be saved in %s.'%os.getcwd()
        outputDir= ''
        
    filename= 'depthMap_%ssigma_%s-band_%sMethod'%(SNRthreshold, band, methodName)
    flatSkyGrid.write_flat_map(outputDir+filename, depth)
    print 'Wrote %s.npz'% filename

    filename= 'depthMap_std_%ssigma_%s-band_%sMethod'%(SNRthreshold, band, methodName)
    flatSkyGrid.write_flat_map(outputDir+filename, depth_std)
    print 'Wrote %s.npz'% filename

#############################################
# code from Javier: /global/projecta/projectdirs/lsst/groups/LSS/DC1/scripts/map_utils.py
# creates 5sigma depth HEALPix maps
def binned_statistic(x, values, func, nbins, range):
    '''The usage is approximately the same as the scipy one
    from https://stackoverflow.com/questions/26783719/effic
    iently-get-indices-of-histogram-bins-in-python'''
    from scipy.sparse import csr_matrix
    r0, r1 = range
    mask = (x > r0) &  (x < r1)
    x = x[mask]
    values = values[mask]
    N = len(values)
    digitized = (float(nbins) / (r1-r0) * (x-r0)).astype(int)
    S = csr_matrix((values, [digitized, np.arange(N)]), shape=(nbins, N))
    return np.array([func(group) for group in np.split(S.data, S.indptr[1:-1])])

def depth_map_snr_nonHP(ra, dec, mags, snr, snrthreshold, flatSkyGrid):
    # not based on healpix, original version modified to use flatmaps
    # also added the functionality to add snr_threshold
    good = np.logical_or(np.logical_not(np.isnan(ra)),np.logical_not(np.isnan(dec)))
    pix_nums = np.array(flatSkyGrid.pos2pix(ra, dec))
    
    map_out = np.zeros(flatSkyGrid.get_size())
    map_var_out = np.zeros(flatSkyGrid.get_size())
    #Binned statistic 2d is awfully slow (because it doesn't use the fact that all bins are equal width
    #median_snr, xed, _, _ = binned_statistic_2d(mags,pix_nums,snr,statistic='median',bins=(50,12*nside**2),
    #                                           range=[(20,30),(0,12*nside**2)])
    #bin_centers = 0.5*xed[1:]+0.5*xed[:-1]
    #depth = bin_centers[np.argmin(np.fabs(median_snr-5),axis=0)]
    
    bin_centers = np.linspace(22+6/30.,28-6/30.,30.)
    for px in np.unique(pix_nums):
        mask = px==pix_nums
        if np.count_nonzero(mask)>0:
            median_snr = binned_statistic(mags[mask],snr[mask],np.nanmedian, nbins=30, range=(22,28))
            std_snr = binned_statistic(mags[mask],snr[mask],np.nanstd, nbins=30, range=(22,28))
            
            mask2 = np.isnan(median_snr)==False
            if np.count_nonzero(mask2)>0:
                depth = bin_centers[mask2][np.argmin(np.fabs(median_snr[mask2]-snrthreshold))]
                std = std_snr[np.argmin(np.fabs(median_snr[mask2]-snrthreshold))]
                map_out[px]=depth
                map_var_out[px]= std
            else:
                map_out[px]=0
                map_var_out[px]= 0
        else:
            map_out[px]=0.
            map_var_out[px]= 0
    return map_out, map_var_out

def desc_method(ra, dec, band, mags, snr, flatSkyGrid, SNRthreshold= 5, plotMaps= True,
                          saveMaps= False, outputDir= None):
    # make a histograms of the S/N in bins of magnitude for all objects in a given pixel
    # define the 5 sigma depth as the magnitude of the histogram whose median S/N is ~5.
    # SNRthreshold= 5 => 5sigma depth. can tweak it.
    SNRthreshold= int(SNRthreshold)
    
    print 'Creating %s-band %ssigma depth maps'%(band, SNRthreshold)
    depth, depth_std= depth_map_snr_nonHP(ra, dec,
                                          mags= mags,
                                          snr= snr,
                                          snrthreshold= SNRthreshold,
                                          flatSkyGrid= flatSkyGrid)    
    # the output has zero in out-of-survey region, so need to mask stuff.
    mask, mask= createMeanStdMaps(ra, dec, quantity= ra,
                                  flatSkyGrid= flatSkyGrid,
                                  plotMaps= False)
    ind= np.where(np.isnan(mask)==True) # out of survey
    depth[ind]= np.nan
    depth_std[ind]= np.nan
    
    if plotMaps:
        flatSkyGrid.view_map(depth, posColorbar= True, 
                             title= "Javi's method\nmean %s-band %sigma depth"%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
        flatSkyGrid.view_map(depth_std, posColorbar= True,
                             title= "Javi's method\nstd %s-band %sigma depth"%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
    
    if saveMaps:
        methodName= 'Javis'
        saveDepthMaps(depth, depth_std, methodName, outputDir, SNRthreshold, band, flatSkyGrid, )
        
    return depth, depth_std

#############################################

def flux_err_method(ra, dec, flux_err, band, flatSkyGrid, SNRthreshold= 5, plotMaps= True,
                    saveMaps= False, outputDir= None):
    # 5sigma Magnitude limit= average of 5*flux_err for all objs in each pixel (and then transformed to magnitude)
    # SNRthreshold= 5 => 5sigma depth.
    
    # Since want mags (mean, std) at the end, need to first run the createMeanStdMaps with 
    # 5flux_error to get the mean flux map which can be converted to fluxea.
    # To get std mags, need to run createMeanStdMaps with 5*flux_error converted to mags and keep only the std.
    
    SNRthreshold= int(SNRthreshold)
    print 'Creating %s-band %ssigma depth maps'%(band, SNRthreshold)
    flux_err= np.array(flux_err)
    
    depth, dontcare = createMeanStdMaps(ra, dec,
                                        quantity= SNRthreshold*flux_err,
                                        flatSkyGrid= flatSkyGrid,
                                        plotMaps= False)
    
    # convert from fluxes to mags
    depth= 10.**(23+6)*depth
    depth[~np.isnan(depth)]= -2.5*np.log10(depth[~np.isnan(depth)])+23.9
    
    # find the std.
    quantity= 10.**(23+6)*SNRthreshold*flux_err
    quanity= -2.5*np.log10(quantity)+23.9
    dontcare, depth_std= createMeanStdMaps(ra, dec,
                                           quantity= quantity,
                                           flatSkyGrid= flatSkyGrid)
    
    if plotMaps:
        flatSkyGrid.view_map(depth, posColorbar= True, 
                             title= 'flux_error method\nmean %s-band %sigma depth'%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
        flatSkyGrid.view_map(depth_std, posColorbar= True,
                             title= 'flux_error method\nstd %s-band %sigma depth'%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
    if saveMaps:
        methodName= 'FluxErr'
        saveDepthMaps(depth, depth_std, methodName, outputDir, SNRthreshold, band, flatSkyGrid)
        
    return depth, depth_std

#############################################
def random_sky_std_method(ra, dec, sky_std, band, flatSkyGrid, SNRthreshold= 5, plotMaps= True,
                          saveMaps= False, outputDir= None):
    # 5sigma depth= -2.5*np.log10(SNRthreshold*sky_std)+27.0,
    # Using the random sky_std as here: https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points-for-dr1/
    # SNRthreshold= 5 => 5sigma depth.
    SNRthreshold= int(SNRthreshold)
    print 'Creating %s-band %ssigma depth maps'%(band, SNRthreshold)
    sky_std= np.array(sky_std)
    depth, depth_std = createMeanStdMaps(ra, dec,
                                         quantity= -2.5*np.log10(SNRthreshold*sky_std)+27.0,
                                         flatSkyGrid= flatSkyGrid,
                                         plotMaps= plotMaps,
                                         quantityName= 'random-sky_std method\n%s-band %ssigma depth'%(band, SNRthreshold))
    if saveMaps:
        methodName= 'randomSkyStd-isPrimary'
        saveDepthMaps(depth, depth_std, methodName, outputDir, SNRthreshold, band, flatSkyGrid, )
        
    return depth, depth_std

#############################################
def depth_map_meanSNRrange(ra, dec, mags, snr, snrthreshold, flatSkyGrid):
    # 5sigma depth= mean of mags of galaxies with 4<SNR<6
    good = np.logical_or(np.logical_not(np.isnan(ra)),np.logical_not(np.isnan(dec)))
    pix_nums = np.array(flatSkyGrid.pos2pix(ra, dec))

    map_out = np.zeros(flatSkyGrid.get_size())
    map_var_out = np.zeros(flatSkyGrid.get_size())
    
    for px in np.unique(pix_nums):
        mask = px==pix_nums
        if np.count_nonzero(mask)>0:
            map_out[px]= np.mean(mags[mask][(snr[mask]>snrthreshold-1) & (snr[mask]<snrthreshold+1)])
            map_var_out[px]= np.std(mags[mask][(snr[mask]>snrthreshold-1) & (snr[mask]<snrthreshold+1)])
        else:
            map_out[px]=0
            map_var_out[px]= 0
            
    return map_out, map_var_out

def dr1paper_method(ra, dec, band, mags, snr, flatSkyGrid, SNRthreshold= 5, plotMaps= True,
                          saveMaps= False, outputDir= None):
    # follow the paper: choose gals with 4<SNR<6 for 5sigma depth.
    # SNRthreshold= 5 => 5sigma depth.
    SNRthreshold= int(SNRthreshold)
    print 'Creating %s-band %ssigma depth maps'%(band, SNRthreshold)
    depth, depth_std= depth_map_meanSNRrange(ra, dec,  mags= mags,
                                             snr= snr,
                                             snrthreshold= SNRthreshold,
                                             flatSkyGrid= flatSkyGrid)
    # the output has zero in out-of-survey region, so need to mask stuff.
    mask, mask= createMeanStdMaps(ra, dec, quantity= ra,
                                  flatSkyGrid= flatSkyGrid,
                                  plotMaps= False)
    ind= np.where(np.isnan(mask)==True) # out of survey
    depth[ind]= np.nan
    depth_std[ind]= np.nan
    
    if plotMaps:
        flatSkyGrid.view_map(depth, posColorbar= True,
                             title= "dr1paper method\nmean %s-band %sigma depth"%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
        flatSkyGrid.view_map(depth_std, posColorbar= True,
                             title= "dr1paper method\nstd %s-band %sigma depth"%(band, SNRthreshold),
                             xlabel='ra', ylabel='dec')
    
    if saveMaps:
        methodName= 'dr1paper'
        saveDepthMaps(depth, depth_std, methodName, outputDir, SNRthreshold, band, flatSkyGrid, )
        
    return depth, depth_std
