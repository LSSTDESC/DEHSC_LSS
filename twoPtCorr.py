import treecorr
import numpy as np
import time

def runTreeCorr(data_ra, data_dec, random_ra, random_dec, # all in degrees
                minSep, maxSep, nBins): 
    print 'Running with %s data pts, %s random pts'%(len(data_ra), len(random_ra))
    startTime= time.time()
    dataCatalog= treecorr.Catalog(ra= data_ra, dec=  data_dec,
                                  ra_units='degrees', dec_units='degrees')
    randCatalog = treecorr.Catalog(ra= random_ra, dec= random_dec,
                                   ra_units='degrees', dec_units='degrees')
    DD = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units='degrees')
    RR= treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units='degrees')
    DR = treecorr.NNCorrelation(min_sep= minSep, max_sep= maxSep, nbins= nBins, sep_units='degrees')

    # use LS
    DD.process(dataCatalog)
    RR.process(randCatalog)
    DR.process(dataCatalog, randCatalog)
    # calculate the correlation function
    wtheta, varxi = DD.calculateXi(RR, DR)
    theta = np.exp(DD.logr) # use bin centers
    wtheta_sig = np.sqrt(varxi)
    
    print 'Time taken: %s (s)'%(time.time()-startTime)

    return theta, wtheta, wtheta_sig
