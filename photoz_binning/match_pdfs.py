# the goal here to match the pdfs from HSC (using the provided object ID) with the galaxies in our catalog
# essentially allows us to retain only the pdfs for the objects in our catalog
import numpy as np
import pandas as pd
import os
from astropy.table import Table
import time
from astropy.io import fits

startTime = time.time()

datapath = '/global/cscratch1/sd/damonge/HSC/HSC_processed' 

##############################################################################################################################
# run over all field
for field in ['aegis', 'gama09h', 'gama15h', 'hectomap', 'vvds', 'wide12h', 'xmmlss']:
    print('\nStarting with %s'%field)
    pdfs_path = '/global/cscratch1/sd/awan/hsc_pdfs/%s'%field

    # ------------------------------------------------------------------------------------------------------------------------
    # read in the cleaned catalog to get the object IDs
    for j, folder in enumerate([f for f in os.listdir(datapath) \
                                if not f.__contains__('.fits') and f.__contains__(field.upper
                                                                                  ())]):
        # read only the Catalog fits file for the field
        for i, file in enumerate([f for f in os.listdir('%s/%s'%(datapath, folder)) if f.__contains__('Catalog')]):
            if i>0 : raise ValueError('Something is wrong. Have more than one Catalog file in %s/%s'%(datapath, folder))

            print('\nReading in %s'%file)
            dat = Table.read('%s/%s/%s'%(datapath, folder,file), format='fits')
        if j==0:
            hscdata = dat.to_pandas()
            print('No concatenation needed. Shape: %s'%(np.shape(hscdata),))
        else:
            raise ValueError('Something is not right. Have more than one folder for %s: %s'%(field, folder))
            #hscdata = pd.concat([hscdata, dat.to_pandas()])
            #print('Data read in. Shape: ', np.shape(dat.to_pandas()))
            #print('Concatenated data. Shape: ', np.shape(hscdata))
        print('')
    
    ##############################################################################################################################
    # read in the pdfs.
    files = [f for f in os.listdir(pdfs_path)]
    print('\n%s patch fits files for %s'%(len(files), field))
    for i, file in enumerate([f for f in files if f.__contains__('.fits')]):  # files for different patches
        print('Reading in %s'%file)
        hdul = fits.open('%s/%s'%(pdfs_path,file)) 
        data_readin = np.array(hdul[1].data)   # pdfs
        bins_readin = np.array(hdul[2].data)   # z_bins
        
        # now need to restructure the data
        if (i==0):  # first patch
            bins = []
            pdfs = {}

        bins_ = []
        for j in range(len(bins_readin)):
            bins_.append(bins_readin[j][0])
        if i==0:
            bins = bins_
        else:
            if bins != bins_:
                raise ValueError('Bins dont match: %s vs. %s'%(bins, bins_))

        for i in range(len(data_readin)):
            pdfs[data_readin[i][0]] = data_readin[i][1]

    bins = np.array(bins)

    ##############################################################################################################################
    # match using IDs
    print('\nMatching IDs now ... ')
    ids, pdfs_cat = [], []

    for i, objID in enumerate(hscdata['object_id']):
        if objID not in pdfs.keys():
            print('\n%s not in the chosen pdfs'%objID)
        else:
            ids.append(objID)
            pdfs_cat.append(pdfs[objID])

    ##############################################################################################################################
    # save data
    hdu = fits.PrimaryHDU(pdfs_cat)
    hdul = fits.HDUList([hdu])
    filename = '/global/cscratch1/sd/awan/hsc_matched_pdfs/matched_pdfs_%s.fits'%field
    hdul.writeto(filename)
    print('\nSaved %s'%filename)

    hdu = fits.PrimaryHDU(ids)
    hdul = fits.HDUList([hdu])
    filename = '/global/cscratch1/sd/awan/hsc_matched_pdfs/matched_ids_%s.fits'%field
    hdul.writeto(filename)
    print('\nSaved %s'%filename)

    hdu = fits.PrimaryHDU(bins)
    hdul = fits.HDUList([hdu])
    filename = '/global/cscratch1/sd/awan/hsc_matched_pdfs/matched_bins_%s.fits'%field
    hdul.writeto(filename)
    print('\nSaved %s'%filename)

    print('\nTime taken (s): ', (time.time()-startTime))