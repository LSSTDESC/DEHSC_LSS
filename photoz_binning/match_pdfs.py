# The goal here to match the pdfs from HSC (using the provided object ID) with the galaxies in our processed catalogs.
#
# Essentially allows us to retain only the pdfs for the objects that remain in our final catalog.
#
##############################################################################################################################
import numpy as np
import pandas as pd
import os
from astropy.table import Table
import time
from astropy.io import fits
import pandas as pd

##############################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--data_main_path', dest='data_main_path',
                  help='Path to the processed data.',
                  default='/global/cscratch1/sd/damonge/HSC/HSC_processed')
parser.add_option('--pdfs_main_path', dest='pdfs_main_path',
                  help='Path to the PDFs.',
                  default='/global/cscratch1/sd/awan/hsc_pdfs/')
parser.add_option('--fields', dest='fields',
                  help='List of fields to consider',
                  default='wide_aegis, wide_gama09h, wide_gama15h, wide_hectomap, wide_vvds, \
                            wide_wide12h, wide_xmmlss, deep_cosmos, deep_elaisn1, deep_xmmlss, deep_deep23')
parser.add_option('--PZalg', dest='PZalg',
                  help='List of PZ algorithms to consider',
                  default='ephor, ephor_ab, demp, frankenz')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where all the matched files should be stored; directory should exist already.',
                  default='/global/cscratch1/sd/awan/hsc_matched_pdfs')

##############################################################################################################################
startTime = time.time()
(options, args) = parser.parse_args()
print('\nOptions: %s'%options)

# read in the inputs
data_main_path = options.data_main_path
pdfs_main_path = options.pdfs_main_path
fields = options.fields
PZalg = options.PZalg
outDir = options.outDir

# format the fields
fields = [f.strip() for f in list(fields.split(','))]
PZalg = [f.strip() for f in list(PZalg.split(','))]
##############################################################################################################################
# run over all field
for field in fields:
    print('\n--------------------------------------------------------------------')
    print('\nStarting with %s'%field)
    # ------------------------------------------------------------------------------------------------------------------------
    # read in the cleaned catalog to get the object IDs
    for j, folder in enumerate([f for f in os.listdir(data_main_path) \
                                if not f.__contains__('.fits') and f.__contains__(field.upper())]):
        # read only the Catalog fits file for the field
        for i, cat_file in enumerate([f for f in os.listdir('%s/%s'%(data_main_path, folder)) if f.__contains__('Catalog')]):
            if i>0 : raise ValueError('Something is wrong. Have more than one Catalog file in %s/%s'%(data_main_path, folder))

            print('\nReading in %s'%cat_file)
            dat = Table.read('%s/%s/%s'%(data_main_path, folder, cat_file), format='fits')
        if j==0:
            hscdata = dat.to_pandas()
            print('Data read in. Shape: %s'%(np.shape(hscdata),))
        else:
            raise ValueError('Something is not right. Have more than one folder for %s: %s'%(field, folder))

        print('')

    # ------------------------------------------------------------------------------------------------------------------------
    for alg in PZalg:
        print('\nStarting with %s'%alg)
        pdfs_path = '%s/%s/%s'%(pdfs_main_path, field, alg)
        ##############################################################################################################################
        # read in the pdfs.
        patch_files = [f for f in os.listdir(pdfs_path) if f.__contains__('.fits')]
        print('\n%s patch fits files for %s'%(len(patch_files), field))
        for i, file in enumerate(patch_files):  # files for different patches
            print('Reading in %s'%file)
            hdul = fits.open('%s/%s'%(pdfs_path, file))
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
        hdr = fits.Header()
        hdr['FIELD'] = field
        hdr['MatchCat'] = cat_file
        hdr['nPatch'] = len(patch_files)
        prm_hdu = fits.PrimaryHDU(header=hdr)
        # data to save
        cat_hdu = fits.table_to_hdu(Table(pd.DataFrame(pdfs_cat).values))
        ids_hdu = fits.table_to_hdu(Table(pd.DataFrame(ids).values))
        bins_hdu = fits.table_to_hdu(Table(pd.DataFrame(bins).values))
        # save it
        hdul = fits.HDUList([prm_hdu, cat_hdu, ids_hdu, bins_hdu])

        filename = '%s/matched_pdfs_ids_bins_%s_%s.fits'%(outDir, field, alg)
        hdul.writeto(filename, overwrite=True)

        print('\nSaved %s'%filename)

        time_taken = time.time()-startTime
        if (time_taken>60.):
            print('\nTime taken %.2f min '%(time_taken/60.) )
        else:
            print('\nTime taken %.2f sec '%time_taken)