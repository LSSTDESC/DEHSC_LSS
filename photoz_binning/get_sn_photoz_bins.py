# The goal here to calculate the SN as a function of Nbin.
#
# We read in the matched PDFs so that we have a PDF for each of the galaxy in the processed catalog. Then, we find bin edges for different Nbin (number of photo-z bins), using the photo-z for each galaxy; these bins are selected s.t. each bin has roughly the same number of galaxies.
#
# Using the point estimates and the PDFs, we calculate dn/dz for each bin; see more in the `get_dn_dz` method. Finally, we use CCL to calculate the signal, the power spectra for each bin pair (i,j), while the point estimates are used to estimate the number density of galaxies in each bin to estimate the shot noise. The final output cell plots the S/N as a function of Nbin.
#
# A few specifics: we work with `AEGIS` field and `ephor_ab` outputs, for $2<\ell<2000$; `pz_mc_eab` are used for photo-z point estimates and Nbin=1-6. Also, for now, the area of the field is assumed to be 108/7 ~ 15.4 deg2.
#
#
##############################################################################################################################
import numpy as np
import pandas as pd
import os
from astropy.table import Table
import time
from astropy.io import fits
from photoz_bin_sn_utils import get_bin_edges, calc_sn
import matplotlib.pyplot as plt
import sys
sys.path.append('%s/..'%os.getcwd())
import flatmaps as fm

##############################################################################################################################
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--cat_data_main_path', dest='data_main_path',
                  help='Path to the processed data.',
                  default='/global/cscratch1/sd/damonge/HSC/HSC_processed')
parser.add_option('--pdfs_main_path', dest='pdfs_main_path',
                  help='Path to the matched PDFs.',
                  default='/global/cscratch1/sd/awan/hsc_matched_pdfs_old/')
parser.add_option('--fields', dest='fields',
                  help='List of fields to consider',
                  default='wide_aegis, wide_gama09h, wide_gama15h, wide_hectomap, wide_vvds, \
                           wide_wide12h, wide_xmmlss, deep_cosmos, deep_elaisn1, deep_xmmlss, deep_deep23')
parser.add_option('--PZalg', dest='PZalg',
                  help='List of PZ algorithms to consider',
                  default='ephor, ephor_ab, demp, frankenz')
parser.add_option('--outDir', dest='outDir',
                  help='Path to the folder where the plots will be saved.',
                  default='/global/cscratch1/sd/awan/lsst_output/hsc_output/')

##############################################################################################################################
startTime = time.time()
(options, args) = parser.parse_args()
print('\nOptions: %s'%options)

# read in the inputs
datapath = options.data_main_path
pdfs_path = options.pdfs_main_path
fields = options.fields
PZalg = options.PZalg
outDir = options.outDir

# format the fields
fields = [f.strip() for f in list(fields.split(','))]
PZalg = [f.strip() for f in list(PZalg.split(','))]

# set up some things
z_phot_key = 'pz_mc_eab'
ell = np.arange(2, 2000)

##############################################################################################################################
# run over all field
for field in fields:
    print('\n--------------------------------------------------------------------')
    print('\nStarting with %s'%field)
    # ------------------------------------------------------------------------------------------------------------------------
    # read in the cleaned catalog.
    for j, folder in enumerate([f for f in os.listdir(datapath) \
                                if not f.__contains__('.fits') and f.__contains__(field.upper())]):
        # read only the Catalog fits file for the field
        for i, cat_file in enumerate([f for f in os.listdir('%s/%s'%(datapath, folder)) if f.__contains__('Catalog')]):
            if i>0 : raise ValueError('Something is wrong. Have more than one Catalog file in %s/%s'%(datapath, folder))

            print('\nReading in %s'%cat_file)
            dat = Table.read('%s/%s/%s'%(datapath, folder, cat_file), format='fits')
        if j==0:
            hscdata = dat.to_pandas()
            print('Data read in. Shape: %s'%(np.shape(hscdata),))
        else:
            raise ValueError('Something is not right. Have more than one folder for %s: %s'%(field, folder))

        ##########################################################################################
        # read the mask fraction file
        for i, maskedfrac_file in enumerate([f for f in os.listdir('%s/%s'%(datapath, folder)) if f.__contains__('MaskedFraction')]):
            if i>0 : raise ValueError('Something is wrong. Have more than one MaskedFraction file in %s/%s'%(datapath, folder))

            print('\nReading in %s'%maskedfrac_file)
            fskb, mskfrac = fm.read_flat_map('%s/%s/%s'%(datapath, folder, maskedfrac_file))
            patch_area = np.sum(mskfrac)*np.radians(fskb.dx)*np.radians(fskb.dy)   # in Sr

    sns, all_bins = {}, {}
    # ------------------------------------------------------------------------------------------------------------------------
    for alg in PZalg:
        print('\nStarting with %s'%alg)
        filename= 'matched_pdfs_ids_bins_wide_%s_%s.fits'%(field, alg)
        print('Reading in %s'%filename)
        hdul = fits.open('%s/%s'%(pdfs_path, filename))

        pdfs = np.array(Table(hdul[1].data))
        pdfs = np.array(pd.DataFrame(pdfs))  # somehow np.array isnt enough to make the table come out nice

        ids = np.array(Table(hdul[2].data)['col0'])
        bins = np.array(Table(hdul[3].data)['col0'])

        # stack the pdfs to estimate N(z)
        n_z = np.sum(pdfs, axis=0)

        # --------------------------------------------------------------
        ### Find the bin edges to consider for different number of bins
        z_phots = []
        n_bin_list = []
        for i in range(1,7):
            n_bin_list.append(i)
            print('----------------------------------\nnbin=%s'%i)
            out = get_bin_edges(nbin=i, hsc_z_phot=hscdata[z_phot_key], z_bins=bins)
            z_phots.append(out)

        # --------------------------------------------------------------
        ### Run CCL and calculate S/N for each Nbin
        sn = np.zeros(len(z_phots))
        for i, z_phot in enumerate(z_phots):
            sn[i] = calc_sn(z_phot=z_phot, z_bins=bins, hsc_z_phot=hscdata[z_phot_key], hsc_ids=hscdata['object_id'],
                            matched_pdf_ids=ids.copy(), matched_pdfs=pdfs.copy(), n_z=n_z, ell=ell, area_in_sr=patch_area,
                            plot_cls=False)

        fsky = patch_area/(4*np.pi)  # total sky area: 4pi Sr
        print('\n## fsky: %s\n'%fsky)
        sns[alg] = (fsky/2.)*np.sqrt(sn)
        all_bins[alg] = z_phots

    # plot SN as a function of Nbin
    plt.clf()
    for key in sns:
        plt.plot(n_bin_list, sns[key], 'o-', label=key)
    plt.xlabel('Nbin')
    plt.ylabel('S/N')
    plt.gca().ticklabel_format(style='sci', scilimits=(-3,4),axis='y')
    plt.xticks(n_bin_list, n_bin_list)#, rotation=70)
    plt.legend()
    plt.show()

    print('S/N: %s'%sns)
    print('all_bins: %s'%all_bins)
    # ------------------------------------------------------------------------------------------------------------------------
    time_taken = time.time()-startTime
    if (time_taken>60.):
        print('\nTime taken %.2f min '%(time_taken/60.) )
    else:
        print('\nTime taken %.2f sec '%time_taken)