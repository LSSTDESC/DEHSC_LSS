import numpy as np
from .map_utils import createMeanStdMaps, createCountsMap

def fluxerr_method(ra, dec, flux_err, fsk, snrthreshold= 5,
                   interpolate=False,count_threshold=4) :
    # 5sigma Magnitude limit= average of 5*flux_err for all objs in each pixel (and then transformed to magnitude)
    # snrthreshold= 5 => 5sigma depth.
    
    # Since want mags (mean, std) at the end, need to first run the createMeanStdMaps with 
    # 5flux_error to get the mean flux map which can be converted to fluxea.
    # To get std mags, need to run createMeanStdMaps with 5*flux_error converted to mags and keep only the std.

    depth, _ = createMeanStdMaps(ra, dec,
                                 quantity= snrthreshold*flux_err,
                                 fsk= fsk)
    
    # convert from fluxes to mags
    depth= 10.**(23+6)*depth
    depth[~np.isnan(depth)]= -2.5*np.log10(depth[~np.isnan(depth)])+23.9
    
    # find the std.
    quantity= 10.**(23+6)*snrthreshold*flux_err
    quanity= -2.5*np.log10(quantity)+23.9
    dontcare, depth_std= createMeanStdMaps(ra, dec,
                                           quantity= quantity,
                                           fsk= fsk)

    #Zeros in empty pixels
    nc=createCountsMap(ra,dec,fsk)
    depth[nc<1]=0
    depth_std[nc<1]=0

    if interpolate:
        from scipy.interpolate import griddata
        idgood=np.where(nc>count_threshold)[0]
        coords_all=np.array(fsk.pix2pos(np.arange(fsk.npix))).T
        depth=griddata(coords_all[idgood],depth[idgood],coords_all,
                       method='nearest',fill_value=0)
        depth_std=griddata(coords_all[idgood],depth_std[idgood],coords_all,
                           method='nearest',fill_value=0)
    return depth, depth_std

#############################################
def dr1_method(ra, dec, mags, snr, fsk, snrthreshold,
               interpolate=False,count_threshold=4):
    good = np.logical_or(np.logical_not(np.isnan(ra)),np.logical_not(np.isnan(dec)))
    pix_nums = np.array(fsk.pos2pix(ra, dec))

    mask=(snr>=snrthreshold-1) & (snr<=snrthreshold+1) & (pix_nums>=0)
    mags_use=mags[mask]
    pix_nums=pix_nums[mask]
    
    n_map = np.bincount(pix_nums,minlength=fsk.npix)
    m_map = np.bincount(pix_nums,weights=mags_use,minlength=fsk.npix)
    m2_map = np.bincount(pix_nums,weights=mags_use**2,minlength=fsk.npix)
    depth=np.zeros(fsk.npix)
    depth_std=np.zeros(fsk.npix)
    pix_good=np.where(n_map>0)[0]
    depth[pix_good]=m_map[pix_good]/n_map[pix_good]
    depth_std[pix_good]=np.sqrt(m2_map[pix_good]/n_map[pix_good]-(m_map[pix_good]/n_map[pix_good])**2)
    if interpolate:
        from scipy.interpolate import griddata
        idgood=np.where(n_map>count_threshold)[0]
        coords_all=np.array(fsk.pix2pos(np.arange(fsk.npix))).T
        depth=griddata(coords_all[idgood],depth[idgood],coords_all,
                       method='nearest',fill_value=0)
        depth_std=griddata(coords_all[idgood],depth_std[idgood],coords_all,
                           method='nearest',fill_value=0)
    return depth, depth_std

def get_depth(method,ra,dec,arr1,arr2,fsk,snrthreshold=5,
                   interpolate=False,count_threshold=4):
    """
    Creates a depth map based on the positions and fluxes of a set of objects.
    :param method: method used to compute the depth map. Allowed values: 'dr1', 'desc' and 'fluxerr'.
    :param ra: right ascension for each object.
    :param dec: declination for each object.
    :param arr1: measurement of the flux (if using 'fluxerr') or magnitude (otherwise) for each object.
    :param arr2: measurement of the S/N for each object (or `None` if using 'fluxerr').
    :param fsk: flatmaps.FlatMapInfo object describing the geometry of the output map.
    :param snrthreshold: S/N cut to use.
    """
    
    snrthreshold= int(snrthreshold)
    print('Creating %ssigma depth maps'%(snrthreshold))
    if method=='dr1' :
        depth,depth_std=dr1_method(ra,dec,mags=arr1,snr=arr2,fsk=fsk,
                                   snrthreshold=snrthreshold,
                                   interpolate=interpolate,
                                   count_threshold=count_threshold)
    elif method=='fluxerr' :
        depth,depth_std=fluxerr_method(ra,dec,flux_err=arr1,fsk=fsk,
                                       snrthreshold=snrthreshold,
                                       interpolate=interpolate,
                                       count_threshold=count_threshold)
    else :
        raise KeyError("Unknown method "+method)

    return depth,depth_std
