import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import copy
import pandas as pd
import seaborn as sns

def getDensity_magCuts(depthmap, band, galDepth, galMags, plotTitle= '', depthRange= [25, 27], deldepthStep= 0.3):
    galDepth= np.array(galDepth)
    galMags= np.array(galMags)
    density= OrderedDict()
    depthBins= np.arange(depthRange[0], depthRange[1]+deldepthStep, deldepthStep)
    
    fontsize= 14
    # plot/store 
    plt.clf()
    fig, axes = plt.subplots(1,2)
    ax= axes[0]
    for i in range(len(depthBins)):
        if (i>0):  # depth bins has left ends
            magCutTag= '%s<=%s<%s'%(depthBins[i-1], band, depthBins[i])
            # find indices for objects with  depth is within the depth range.
            magCutInd= np.where((galDepth>= depthBins[i-1]) & (galDepth < depthBins[i]))[0]
            
            # weight by the number of pixels.
            weight= np.ones(len(magCutInd))/np.count_nonzero((depthmap>=depthBins[i-1])&(depthmap<depthBins[i]))
            #  histrogram the objects' mags that satisfy the 5sigma depth cut.
            density[magCutTag]= ax.hist(galMags[magCutInd], weights= weight,
                                         label= magCutTag, histtype= 'step', alpha= 1., lw= 3, bins= 50)
    ax.legend(loc= 'upper left')
    ax.set_xlabel('mag', fontsize= fontsize)
    ax.set_ylabel('Number of objects/Number of pixels', fontsize= fontsize)
    #ax.set_xlim(0, maxDepth+5*step)
    #plt.gcf().set_size_inches(10,6)
    #plt.show()
    
    # plot cumulative sums
    # try to look at the plot above another way?
    ax= axes[1]
    for cut in density.keys():    
        ax.plot(density[cut][1][1:], np.cumsum(density[cut][0]), label= cut)
    ax.legend(loc= 'upper left')
    ax.set_xlabel('Mag', fontsize= fontsize)
    ax.set_ylabel('CumSum(Number of objects/Number of pixels)', fontsize= fontsize)
    ax.set_xlim(22,30)
    
    for ax in axes:
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)

    plt.suptitle(plotTitle, fontsize= fontsize, fontweight="bold",)
    plt.xlim(22,30)
    fig.set_size_inches(20,6)
    plt.xlim(22,30)
    plt.show()

    return density
    

def getMask(depth_map, band, limitingMag, plotMap= True, flatSkyGrid= None, title= ''):
    print 'Finding mask for %s<%s'%(band, limitingMag)
    depth_map= np.array(depth_map)
    ind= np.where(depth_map> limitingMag)[0]
    mask= copy.deepcopy(depth_map)
    
    mask[~np.isnan(mask)]= 0.  # keep the mask
    mask[ind]= 1.
    
    if plotMap:
        if flatSkyGrid is None:
            raise ValueError('Need flatSkyGrid to plot here.')
        flatSkyGrid.view_map(mask, posColorbar= False, 
                             title= "%s\n mask: %s<%s"%(title,  band, limitingMag),
                             xlabel='ra', ylabel='dec')
    return mask


def getRandomCatalog(flatSkyGrid, surveyMask, minRA, maxRA, minDec, maxDec, nData,  plotMap= True):
    # ra, dec in degrees
    # uniform random catalog. no geometry yet.
    rand= pd.DataFrame()
    rand['ra'] = np.rad2deg(np.random.uniform(np.deg2rad(minRA), np.deg2rad(maxRA), nData*2))
    rand['sinDec']= np.random.uniform(np.deg2rad(minDec), np.deg2rad(maxDec), nData*2)
    rand['dec'] = np.rad2deg(np.arcsin(rand['sinDec']))
 
    # create random catalog with geometry: find pixels coresponding to the random ra, dec. keep only those that
    # have depth= 1 in the depth map.
    # find pixel numbers corresponding to the ra, decs.
    pixelNums= flatSkyGrid.pos2pix(rand['ra'], rand['dec'])
    # find pixel numbers in survey
    pixelNumsInSurvey= np.where(surveyMask>0)[0]

    # find out whether the random-pixelNum is in survey. mask= True <=> yes
    mask= np.in1d(pixelNums, pixelNumsInSurvey)
    badind= np.where(mask==False)[0]

    randCatalog= rand.drop(badind, axis= 0)
    print '%s entries in the random catalog.'%len(randCatalog)

    if plotMap:
        plt.clf()
        sns.jointplot(x=randCatalog['ra'], y=randCatalog['dec'], kind="hex", color="k")
        plt.show()
        
    return randCatalog
