import numpy as np
from scipy.ndimage import label

def createCountsMap(ra, dec, flatSkyGrid, returnMap= True, plotMap= False, quantityName= ''):
    flatmap= flatSkyGrid.pos2pix(ra, dec)
    mp= np.bincount(flatmap, weights= None, minlength= flatSkyGrid.get_size())

    if plotMap:
        flatSkyGrid.view_map(mp, posColorbar= True, title= '%s count'%quantityName, xlabel='ra', ylabel='dec')
    if returnMap:
        return mp

def createMeanStdMaps(ra, dec, quantity, flatSkyGrid, returnMaps= True, plotMaps= False, quantityName= '',nan_outside=False):
    pix_ids=flatSkyGrid.pos2pix(ra, dec)
    mp= np.bincount(pix_ids, weights= None, minlength= flatSkyGrid.get_size())
    mpWeighted= np.bincount(pix_ids, weights= quantity, minlength= flatSkyGrid.get_size())
    mpWeightedSq= np.bincount(pix_ids, weights= quantity**2, minlength= flatSkyGrid.get_size())
    idgood=np.where(mp>0)[0];
    mean=np.zeros(len(mp)); std=np.zeros(len(mp))
    mean[idgood]= mpWeighted[idgood]/mp[idgood]
    #std[idgood]= np.sqrt(np.fabs((mpWeightedSq[idgood]/mp[idgood])-mean[idgood]**2))/np.sqrt(mp[idgood]+0.)
    std[idgood]= np.sqrt(np.fabs(((mpWeightedSq[idgood]/mp[idgood])-mean[idgood]**2)/(mp[idgood]+0.)))
    idbad=np.where(mp<=0)[0]
    if nan_outside :
        mean[idbad]=np.nan
        std[idbad]=np.nan
    else :
        mean[idbad]=0
        std[idbad]=0

    if plotMaps:
        flatSkyGrid.view_map(mean, posColorbar= True, title= 'mean %s'%quantityName, xlabel='ra', ylabel='dec')
        flatSkyGrid.view_map(std, posColorbar= True, title= 'std %s'%quantityName, xlabel='ra', ylabel='dec')
    if returnMaps:
        return mean, std

def removeDisconnected(mp,fsk) :
    labeled_array,num_features=label(mp.reshape([fsk.ny,fsk.nx]))
    i0=np.argmax(np.histogram(labeled_array,bins=num_features+1,range=[-0.5,num_features+0.5])[0][1:])+1
    mask_clean=labeled_array.copy().flatten()
    mpo=mp.copy()
    mpo[mask_clean!=i0]=0

    return mpo
