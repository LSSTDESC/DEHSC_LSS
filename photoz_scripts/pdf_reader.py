from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits
from tqdm import tqdm
import os
from scipy.optimize import minimize
from astropy.table import Table

class hsc_reader:

    """
        UPDATE
    """

    def __init__(self, inputfp = './data/deep_cats/DEEP_COSMOS_Catalog_i24.50.fits'):

        catlist = sorted(glob(inputfp))

        self.object_id = []
        self.ra = []
        self.dec = []

        self.mag_g = []
        self.mag_r = []
        self.mag_i = []
        self.mag_z = []
        self.mag_y = []

        # self.demp_mc = []
        # self.ephor_mc = []
        self.ephor_ab_mc = []
        self.frankenz_mc = []
        self.nnpz_mc = []

        self.catname = []

        for x in tqdm(xrange(len(catlist))):

            this_data = fits.open(catlist[x])[1].data
            self.object_id = self.object_id + list(this_data['object_id'])
            self.ra = self.ra + list(this_data['ra'])
            self.dec = self.dec + list(this_data['dec'])

            self.mag_g = self.mag_g + list(this_data['gcmodel_flux'])
            self.mag_r = self.mag_r + list(this_data['rcmodel_flux'])
            self.mag_i = self.mag_i + list(this_data['icmodel_flux'])
            self.mag_z = self.mag_z + list(this_data['zcmodel_flux'])
            self.mag_y = self.mag_y + list(this_data['ycmodel_flux'])

            # self.demp_mc = self.demp_mc + list(this_data['pz_mc_dem'])
            # self.ephor_mc = self.ephor_mc + list(this_data['pz_mc_eph'])
            self.ephor_ab_mc = self.ephor_ab_mc + list(this_data['pz_mc_eab'])
            self.frankenz_mc = self.frankenz_mc + list(this_data['pz_mc_frz'])
            self.nnpz_mc = self.nnpz_mc + list(this_data['pz_mc_nnz'])

            self.catname = self.catname + [catlist[x].split('_')[-3],]*len(this_data['object_id'])

        self.object_id = np.array(self.object_id, dtype = int)
        self.ra = np.array(self.ra, dtype = float)
        self.dec = np.array(self.dec, dtype = float)

        self.mag_g = np.array(self.mag_g, dtype = float)
        self.mag_r = np.array(self.mag_r, dtype = float)
        self.mag_i = np.array(self.mag_i, dtype = float)
        self.mag_z = np.array(self.mag_z, dtype = float)
        self.mag_y = np.array(self.mag_y, dtype = float)

        # self.demp_mc = np.array(self.demp_mc, dtype = float)
        # self.ephor_mc = np.array(self.ephor_mc, dtype = float)
        self.ephor_ab_mc = np.array(self.ephor_ab_mc, dtype = float)
        self.frankenz_mc = np.array(self.frankenz_mc, dtype = float)
        self.nnpz_mc = np.array(self.nnpz_mc, dtype = float)

        self.catname = np.array(self.catname, dtype = str)



def fits_converter(inputfp, outputfp):

    hdulist = fits.open(inputfp)
    pdfs = hdulist[1].data
    pdf_cols = hdulist[1].data.names
    ids = hdulist[2].data
    bins = hdulist[3].data

    # Restructure the pdfs as an array

    new_pdfs = np.vstack([pdfs[pdf_cols[x]] for x in xrange(len(pdf_cols))]).T 

    # Stuff everything into Fits data types

    hdr = fits.Header()
    primary_hdu = fits.PrimaryHDU(header=hdr)
    
    # The IDs and pdfs can go in the same table because they are both length N

    col1 = fits.Column(name = 'object_id', format = 'K', array = np.array(ids, dtype = int))
    col2 = fits.Column(name = 'pdf', format = '601E', array = new_pdfs)
    cols = fits.ColDefs([col1, col2])
    pdf_hdu = fits.BinTableHDU.from_columns(cols)

    # The bins go in a different table because we only need to store them once

    bincol = fits.Column(name = 'bins', format = 'E', array = np.array(bins, dtype = float))
    bincols = fits.ColDefs([bincol])
    bin_hdu = fits.BinTableHDU.from_columns(bincols)

    hdulist = fits.HDUList([primary_hdu, pdf_hdu, bin_hdu])

    hdulist.writeto(outputfp)




class reader:

    # def __init__(self, name = 'ephor', cachefolder = './data/cache/'):

    #     self.name = name

    #     inputfp = './data/pdr1_' + name + '_deep_cosmos/'

    #     file_list = sorted(glob(inputfp + '*'))

    #     self.object_id = np.array([])
    #     self.pdf = []
    #     self.bins = np.array([])

    #     idlist = []
    #     pdflist = []
    #     binlist = []
    #     hsc_indices = []

    #     # This loop goes through each of the files in the PDF directory

    #     for filename in tqdm(file_list):

    #         temp_hsc_indices = []
    #         temp_match_indices = []

    #         hdulist = fits.open(filename)

    #         # Set the cache file name for every file in the PDF directory

    #         cachefile = cachefolder + self.name + '_' + filename.split('/')[-1]
    #         cachefile = cachefile[:-5] + '.dat'

    #         # If the cache file doesn't exist, create it

    #         if not os.path.isfile(cachefile):
    #             for x in tqdm(xrange(len(hsc_cat.object_id))):

    #                 # Find which IDs in the PDF file match HSC IDs

    #                 index = np.where(hsc_cat.object_id[x] == hdulist[1].data['ID'])[0]

    #                 # If there is an ID match, add it to hsc_indices and match_indices

    #                 if len(index) != 0:
    #                     index = index[0]
    #                     temp_hsc_indices.append(x)
    #                     temp_match_indices.append(index)

    #             tabledata = np.array((temp_hsc_indices, temp_match_indices)).T
    #             if not os.path.isdir('./data/cache/'):
    #                 os.mkdir('./data/cache/')
    #             np.savetxt(cachefile, tabledata, header = 'hsc_index pdf_index \n', fmt = '  %8i')

    #         else:
    #             temp_hsc_indices, temp_match_indices = np.loadtxt(cachefile, unpack = True, dtype = int)

    #         hsc_indices = hsc_indices + list(temp_hsc_indices)

    #         # Extract only the PDFs which have matching IDs with HSC

    #         pdflist.append(hdulist[1].data['PDF'][temp_match_indices])
    #         if self.name == 'frankenz':
    #             binlist.append(hdulist[2].data['zgrid'])
    #         else:
    #             binlist.append(hdulist[2].data['BINS'])

    #         if all([all(binlist[0] == rest) for rest in binlist]):
    #             self.bins = binlist[0]

    #     # Store an attribute which says which objects in the PDF catalog correspond with what objects in the HSC catalog

    #     self.hsc_indices = np.array(hsc_indices)

    #     # Turn self.pdf in to an array composed of many PDFs along axis 0

    #     self.pdf = np.vstack(pdflist)



    def __init__(self, inputfp = './data/match_pdfs/*cosmos*'):

        if inputfp != './data/match_pdfs/*cosmos*':
            print 'WARNING: This code may have bugs when not run on one catalog.'

        self.hsc_cat = hsc_reader()

        file_list = sorted(glob(inputfp))

        self.demp_pdf = []
        self.ephor_pdf = []
        self.ephor_ab_pdf = []
        self.frankenz_pdf = []

        self.demp_bins = None
        self.ephor_bins = None
        self.ephor_ab_bins = None
        self.frankenz_bins = None

        for x in tqdm(xrange(len(file_list))):

            hdulist = fits.open(file_list[x])
            pdfs = hdulist[1].data
            ids = hdulist[2].data
            bins = hdulist[3].data

            if file_list[x].split('_')[-1] == 'demp.fits':
                # Add to demp pdfs
                # self.demp_pdf = self.demp_pdf + list(pdfs)
                self.demp_pdf = pdfs
                if self.demp_bins == None:
                    self.demp_bins = bins

            elif file_list[x].split('_')[-1] == 'ephor.fits':
                # Add to ephor pdfs
                # self.ephor_pdf = self.ephor_pdf + list(pdfs)
                self.ephor_pdf = pdfs
                if self.ephor_bins == None:
                    self.ephor_bins = bins

            elif file_list[x].split('_')[-1] == 'ab.fits':
                # Add to ephor_ab pdfs
                # self.ephor_ab_pdf = self.ephor_ab_pdf + list(pdfs)
                self.ephor_ab_pdf = pdfs
                if self.ephor_ab_bins == None:
                    self.ephor_ab_bins = bins

            elif file_list[x].split('_')[-1] == 'frankenz.fits':
                # Add to frankenz pdfs
                # self.frankenz_pdf = self.frankenz_pdf + list(pdfs)
                self.frankenz_pdf = pdfs
                if self.frankenz_bins == None:
                    self.frankenz_bins = bins

        # self.demp_pdf = np.array(self.demp_pdf)
        # self.ephor_pdf = np.array(self.ephor_pdf)
        # self.ephor_ab_pdf = np.array(self.ephor_ab_pdf)
        # self.frankenz_pdf = np.array(self.frankenz_pdf)



    def get_integrals(self, codename = 'demp', stretch = 1., shift = 0.):

        if codename == 'ephor':
            z_mc = self.hsc_cat.ephor_mc
            pdf = self.ephor_pdf
            bins = self.ephor_bins
        elif codename == 'ephor_ab':
            z_mc = self.hsc_cat.ephor_ab_mc
            pdf = self.ephor_ab_pdf
            bins = self.ephor_ab_bins
        elif codename == 'demp':
            z_mc = self.hsc_cat.demp_mc
            pdf = self.demp_pdf
            bins = self.demp_bins
        elif codename == 'frankenz':
            z_mc = self.hsc_cat.frankenz_mc
            pdf = self.frankenz_pdf
            bins = self.frankenz_bins
        # elif codename == 'mizuki':
        #     z_mc = hsc_cat.mizuki_mc
        # elif codename == 'mlz':
        #     z_mc = hsc_cat.mlz_mc
        # elif codename == 'nnpz':
        #     z_mc = hsc_cat.nnpz_mc

        integrals = []

        # Calculate the cumulative integral from the lower bound to the estimated redshift

        for x in tqdm(xrange(len(pdf))):

            mod_pdf = np.array(pdf[x])
            mod_bins = bins
            this_montecarlo = z_mc[x]

            if stretch != 1.:
                mod_bins = ((mod_bins - this_montecarlo) * stretch) + this_montecarlo
                mod_pdf = mod_pdf[mod_bins > 0]
                mod_bins = mod_bins[mod_bins > 0]

            if shift != 0.:
                mod_bins = mod_bins + shift

            part = np.trapz(mod_pdf[mod_bins < this_montecarlo], x = mod_bins[mod_bins < this_montecarlo])
            total = np.trapz(mod_pdf, x = mod_bins)
            integrals.append(part/total)

        integrals = np.array(integrals)

        return integrals



    def plot_PIT(self, codename, stretch = 1., shift = 0.):

        integrals = self.get_integrals(codename, stretch, shift)

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.hist(integrals, bins = 50, normed = True, range = [0,1], histtype = 'step', color = 'k')
        sp.plot([0,1],[1,1], color = 'r')
        sp.set_xlim(0,1)
        sp.set_ylim(0,5)



    def plot_PIT_fit(self, codename):

        fig = plt.figure(figsize = (16,8))
        sp1 = fig.add_subplot(121)
        sp2 = fig.add_subplot(122)

        integrals = self.get_integrals(codename)

        heights, bins = np.histogram(integrals, bins = 50, range = [0,1], normed = True)

        sp1.step(bins[:-1] + 0.5*(bins[1]-bins[0]), heights, color = 'k')
        sp1.plot([0,1],[1,1], color = 'r')
        sp1.set_xlim(0,1)
        sp1.set_ylim(0,5)


        def calchist((stretch, shift)):

            mod_integrals = self.get_integrals(codename, stretch, shift)
            mod_heights, mod_bins = np.histogram(mod_integrals, bins = 50, range = [0,1], normed = True)

            return sum((mod_heights-1.)**2.)

        result = minimize(calchist, [1,0], method = 'Nelder-Mead')

        fit_stretch, fit_shift = result.x

        print result

        new_integrals = self.get_integrals(codename, fit_stretch, fit_shift)

        new_heights, new_bins = np.histogram(new_integrals, bins = 50, range = [0,1], normed = True)

        sp2.step(new_bins[:-1] + 0.5*(new_bins[1] - new_bins[0]), new_heights, color = 'k')
        sp2.plot([0,1],[1,1], color = 'r')
        sp2.set_xlim(0,1)
        sp2.set_ylim(0,5)

        sp1.set_ylabel('Norm Frequency', family = 'Roboto', weight = 'light', fontsize = 20)

        sp2.text(0.02, 0.98, 'Stretch: %.3f\nShift: %.3f' % (fit_stretch, fit_shift), family = 'Roboto', weight = 'light', fontsize = 20, ha = 'left', va = 'top', transform = sp2.transAxes)
        fig.text(0.5,0.98, self.name.capitalize(), ha = 'center', va = 'center', family = 'Roboto', weight = 'light', fontsize = 24)
        fig.text(0.5, 0.02, 'PIT', family = 'Roboto', weight = 'light', fontsize = 20)

        sp2.set_yticklabels([])

        plt.subplots_adjust(wspace = 0)


