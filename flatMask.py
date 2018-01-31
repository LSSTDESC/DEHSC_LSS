import numpy as np
from scipy.ndimage import label
from createMaps import createCountsMap

def createMask(ra,dec,flags,flatsky_base,reso_mask) :
    """
    Creates a mask based on the position of random objects and a set of flags.
    ra,dec : angular coordinates of the random objects
    flags : list of arrays containing the flags used to mask areas of the sky.
            pixels containing objects with any flags=True will be masked.
            Pass [] if you just want to define a mask based on object positions.
    flatsky_base : FlatMapInfo for the base mask, defined by the presence of not
                   of object in pixels defined by this FlatMapInfo
    reso_mask : resolution of the final mask (dx or dy)
    returns mask and associated FlatMapInfo
    """
    fsg0=flatsky_base

    #Create mask based on object positions
    mpr=createCountsMap(ra,dec,fsg0,returnMap=True,plotMap=False)
    mskr=np.zeros(fsg0.get_size()); mskr[mpr>0]=1

    if(np.sum(mpr*mskr)/np.sum(mskr)<5) :
        raise Warning('Base resolution may be too high %.1lf'%(np.sum(mpr*mskr)/np.sum(mskr)))
    
    if reso_mask>fsg0.dx :
        fsg,mpn=fsg0.d_grade(mpr,int(reso_mask/fsg0.dx+0.5))
    else :
        fsg,mpn=fsg0.u_grade(mpr,int(fsg0.dx/reso_mask+0.5))

    mskn=np.zeros(fsg.get_size()); mskn[mpn>0]=1
    ipix=fsg.pos2pix(ra,dec)
    for flag in flags :
        ipixmask=ipix[flag]
        p=np.unique(ipixmask)
        mskn[p]=0

    #Classify all connected regions
    msk2d=mskn.astype(int).reshape([fsg.ny,fsg.nx])
    labeled_array,num_features=label(msk2d)
    #Identify largest connected region (0 is background)
    i0=np.argmax(np.histogram(labeled_array,bins=num_features+1,range=[-0.5,num_features+0.5])[0][1:])+1
    #Remove all other regions
    mask_clean=labeled_array.copy().flatten()
    mask_clean[mask_clean!=i0]=0
    msk_out=mask_clean.astype(float)
    
    return msk_out,fsg
