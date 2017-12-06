import numpy as np

def createCountsMap(ra, dec, flatSkyGrid, returnMap= True, plotMap= False, quantityName= ''):
    flatmap= flatSkyGrid.pos2pix(ra, dec)
    mp= np.bincount(flatmap, weights= None, minlength= flatSkyGrid.get_size())

    if plotMap:
        flatSkyGrid.view_map(mp, posColorbar= True, title= '%s count'%quantityName, xlabel='ra', ylabel='dec')
    if returnMap:
        return mp

def createMeanStdMaps(ra, dec, quantity, flatSkyGrid, returnMaps= True, plotMaps= False, quantityName= ''):
    flatmap=flatSkyGrid.pos2pix(ra, dec)
    mp= np.bincount(flatmap, weights= None, minlength= flatSkyGrid.get_size())
    mpWeighted= np.bincount(flatmap, weights= quantity, minlength= flatSkyGrid.get_size())
    mpWeightedSq= np.bincount(flatmap, weights= quantity**2, minlength= flatSkyGrid.get_size())
    mean= mpWeighted/mp
    std= np.sqrt((mpWeightedSq/mp)-mean**2)/np.sqrt(mp)
    if plotMaps:
        flatSkyGrid.view_map(mean, posColorbar= True, title= 'mean %s'%quantityName, xlabel='ra', ylabel='dec')
        flatSkyGrid.view_map(std, posColorbar= True, title= 'std %s'%quantityName, xlabel='ra', ylabel='dec')
    if returnMaps:
        return mean, std
