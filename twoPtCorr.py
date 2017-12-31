from __future__  import print_function
import treecorr
import numpy as np
import time
import multiprocessing
def runTreeCorr(data_ra, data_dec, random_ra, random_dec, # all in degrees
                minSep, maxSep, nBins): 
    print('Running with %s data pts, %s random pts'%(len(data_ra), len(random_ra)))
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
    
    print('Time taken: %s (s)'%(time.time()-startTime))

    return theta, wtheta, wtheta_sig

def get_labels(data_ra,data_dec,rnd_ra,rnd_dec,ncen):
    import kmeans_radec
    X= np.zeros((len(data_ra),2))
    X[:,0] = data_ra
    X[:,1] = data_dec
    X2 = np.zeros((len(rnd_ra),2))
    X2[:,0] = rnd_ra
    X2[:,1] = rnd_dec
    cen_guess = np.zeros((ncen,2))
    xmin = np.nanmin(data_ra)
    xmax = np.nanmax(data_ra)
    ymin = np.nanmin(data_dec)
    ymax = np.nanmin(data_dec)
    cen_guess[:,0] = np.linspace(xmin,xmax,ncen)
    cen_guess[:,1] = np.linspace(ymin,ymax,ncen)
    km = kmeans_radec.KMeans(cen_guess,tol=1e-3)
    km.run(X, maxiter=100)
    labels = km.find_nearest(X)
    labels_rnd = km.find_nearest(X2)
    return labels, labels_rnd
def get_labels_one(data_ra,data_dec,ncen):
    import kmeans_radec
    X= np.zeros((len(data_ra),2))
    X[:,0] = data_ra
    X[:,1] = data_dec
    cen_guess = np.zeros((ncen,2))
    xmin = np.nanmin(data_ra)
    xmax = np.nanmax(data_ra)
    ymin = np.nanmin(data_dec)
    ymax = np.nanmin(data_dec)
    cen_guess[:,0] = np.linspace(xmin,xmax,ncen)
    cen_guess[:,1] = np.linspace(ymin,ymax,ncen)
    km = kmeans_radec.KMeans(cen_guess,tol=1e-3)
    km.run(X, maxiter=100)
    labels = km.find_nearest(X) 
    return labels
