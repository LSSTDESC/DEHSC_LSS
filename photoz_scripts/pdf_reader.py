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

frankenz_PIT_fit = [0.853, -0.009]
ephor_ab_PIT_fit = [0.682, -0.016]

class hsc_reader:

    """
        UPDATE
    """

    def __init__(self, inputfp = './data/deep_cats/DEEP_COSMOS_Catalog_i24.50.fits', specfp = './data/11473.csv'):

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

        self.ephor_ab_peak = []
        self.frankenz_peak = []
        self.nnpz_peak = []

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

            self.ephor_ab_peak = self.ephor_ab_peak + list(this_data['pz_mode_eab'])
            self.frankenz_peak = self.frankenz_peak + list(this_data['pz_mode_frz'])
            self.nnpz_peak = self.nnpz_peak + list(this_data['pz_mode_nnz'])

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

        self.ephor_ab_peak = np.array(self.ephor_ab_peak, dtype = float)
        self.frankenz_peak = np.array(self.frankenz_peak, dtype = float)
        self.nnpz_peak = np.array(self.nnpz_peak, dtype = float)

        self.catname = np.array(self.catname, dtype = str)

        self.z_spec = np.ones(len(self.object_id))*-99

        spec_id, spec_z = pd.read_csv(specfp, usecols = ['# object_id', 'specz_redshift']).values.T

        if not os.path.isdir('./data/hsc_cache/'):
            os.mkdir('./data/hsc_cache/')
        if os.path.isfile('./data/hsc_cache/specz.txt'):
            self.z_spec = np.loadtxt('./data/hsc_cache/specz.txt', unpack = True)
        else:

            for x in tqdm(xrange(len(spec_id))):

                self.z_spec[self.object_id == spec_id[x]] = spec_z[x]

            np.savetxt('./data/hsc_cache/specz.txt', self.z_spec.T, fmt = '%.4e')



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
    col2 = fits.Column(name = 'pdf', format = '%iE' % len(bins), array = new_pdfs)
    cols = fits.ColDefs([col1, col2])
    pdf_hdu = fits.BinTableHDU.from_columns(cols)

    # The bins go in a different table because we only need to store them once

    bincol = fits.Column(name = 'bins', format = 'E', array = np.array(bins, dtype = float))
    bincols = fits.ColDefs([bincol])
    bin_hdu = fits.BinTableHDU.from_columns(bincols)

    hdulist = fits.HDUList([primary_hdu, pdf_hdu, bin_hdu])

    hdulist.writeto(outputfp)




class reader:


    def __init__(self, inputfp = './data/readable_match_pdfs/*cosmos*'):

        if inputfp != './data/readable_match_pdfs/*cosmos*':
            print 'WARNING: This code may have bugs when not run on one catalog.'

        self.hsc_cat = hsc_reader()

        file_list = sorted(glob(inputfp))

        self.demp_pdf = []
        self.ephor_pdf = []
        self.ephor_ab_pdf = []
        self.frankenz_pdf = []
        self.nnpz_pdf = []

        self.demp_bins = None
        self.ephor_bins = None
        self.ephor_ab_bins = None
        self.frankenz_bins = None
        self.nnpz_bins = None

        for x in tqdm(xrange(len(file_list))):

            hdulist = fits.open(file_list[x])
            pdfs = hdulist[1].data['pdf']
            ids = hdulist[1].data['object_id']
            bins = hdulist[2].data['bins']

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

            elif file_list[x].split('_')[-1] == 'nnpz.fits':
                # Add to nnpz pdfs
                # self.nnpz_pdf = self.nnpz_pdf + list(pdfs)
                self.nnpz_pdf = pdfs
                if self.nnpz_bins == None:
                    self.nnpz_bins = bins

        # self.demp_pdf = np.array(self.demp_pdf)
        # self.ephor_pdf = np.array(self.ephor_pdf)
        # self.ephor_ab_pdf = np.array(self.ephor_ab_pdf)
        # self.frankenz_pdf = np.array(self.frankenz_pdf)



    def get_integrals(self, codename = 'demp', stretch = 1., shift = 0., ztype = 'zspec', widelims = True):


        if codename == 'ephor_ab':
            if ztype == 'zmc':
                redshift = self.hsc_cat.ephor_ab_mc
            elif ztype == 'zspec':
                redshift = self.hsc_cat.z_spec
            elif ztype == 'zpeak':
                redshift = self.hsc_cat.ephor_ab_peak #Insert peak redshift
            pdf = self.ephor_ab_pdf
            bins = self.ephor_ab_bins
        elif codename == 'frankenz':
            if ztype == 'zmc':
                redshift = self.hsc_cat.frankenz_mc
            elif ztype == 'zspec':
                redshift = self.hsc_cat.z_spec
            elif ztype == 'zpeak':
                redshift = self.hsc_cat.frankenz_peak #Insert peak redshift
            pdf = self.frankenz_pdf
            bins = self.frankenz_bins

        if widelims:
            # Implement wide field magnitude limits
            lim_indices = np.where((self.hsc_cat.mag_r < 26) & (redshift > 0))
        else:
            lim_indices = np.where(redshift > 0)

        redshift = redshift[lim_indices]
        pdf = pdf[lim_indices]

        integrals = []

        # Calculate the cumulative integral from the lower bound to the estimated redshift

        for x in tqdm(xrange(len(pdf))):

            mod_pdf = np.array(pdf[x])
            mod_bins = bins
            this_montecarlo = redshift[x]

            if stretch != 1.:
                mod_bins = ((mod_bins - this_montecarlo) * stretch) + this_montecarlo
                mod_pdf = mod_pdf[mod_bins > 0]
                mod_bins = mod_bins[mod_bins > 0]

            if shift != 0.:
                mod_bins = mod_bins + shift
                mod_pdf = mod_pdf[mod_bins > 0]
                mod_bins = mod_bins[mod_bins > 0]

            part = np.trapz(mod_pdf[mod_bins < this_montecarlo], x = mod_bins[mod_bins < this_montecarlo])
            total = np.trapz(mod_pdf, x = mod_bins)
            integrals.append(part/total)

        integrals = np.array(integrals)

        return integrals


    def get_integrals2(self, codename = 'demp', power = 1., bias = 0., shift = 0., ztype = 'zspec', widelims = True):

        # This version calculates integrals using power multiplication, bias, and shift

        if codename == 'ephor_ab':
            if ztype == 'zmc':
                redshift = self.hsc_cat.ephor_ab_mc
            elif ztype == 'zspec':
                redshift = self.hsc_cat.z_spec
            elif ztype == 'zpeak':
                redshift = self.hsc_cat.ephor_ab_peak #Insert peak redshift
            pdf = self.ephor_ab_pdf
            bins = self.ephor_ab_bins
        elif codename == 'frankenz':
            if ztype == 'zmc':
                redshift = self.hsc_cat.frankenz_mc
            elif ztype == 'zspec':
                redshift = self.hsc_cat.z_spec
            elif ztype == 'zpeak':
                redshift = self.hsc_cat.frankenz_peak #Insert peak redshift
            pdf = self.frankenz_pdf
            bins = self.frankenz_bins

        if widelims:
            # Implement wide field magnitude limits
            lim_indices = np.where((self.hsc_cat.mag_r < 26) & (redshift > 0))
        else:
            lim_indices = np.where(redshift > 0)

        redshift = redshift[lim_indices]
        pdf = pdf[lim_indices]

        integrals = []

        # Calculate the cumulative integral from the lower bound to the estimated redshift

        for x in tqdm(xrange(len(pdf))):

            mod_pdf = np.array(pdf[x])
            mod_bins = bins
            thisredshift = redshift[x]

            if power != 1.:
                mod_pdf = np.power(mod_pdf, power)

            if bias != 0.:
                multipliers = bias * (mod_bins - thisredshift) + 1.
                multipliers[multipliers < 0] = 0
                mod_pdf = mod_pdf * multipliers

            if shift != 0.:
                mod_bins = mod_bins + shift
                mod_pdf = mod_pdf[mod_bins > 0]
                mod_bins = mod_bins[mod_bins > 0]


            part = np.trapz(mod_pdf[mod_bins < thisredshift], x = mod_bins[mod_bins < thisredshift])
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



    def plot_PIT_fit(self, codename, new_integration = True):

        integrals = self.get_integrals(codename)

        heights, bins = np.histogram(integrals, bins = 50, range = [0,1], normed = True)

        if not new_integration:

            def calchist((stretch, shift)):

                mod_integrals = self.get_integrals(codename, stretch, shift)
                mod_heights, mod_bins = np.histogram(mod_integrals, bins = 50, range = [0,1], normed = True)

                return sum((mod_heights-1.)**2.)

            result = minimize(calchist, [1,0], method = 'Nelder-Mead')

        else:

            def calchist((power, bias, shift)):

                mod_integrals = self.get_integrals2(codename, power, bias, shift)
                mod_heights, mod_bins = np.histogram(mod_integrals, bins = 50, range = [0,1], normed = True)

                return sum((mod_heights-1.)**2.)

            result = minimize(calchist, [1,0,0], method = 'Nelder-Mead')

        fit_params = tuple(result.x)

        print result


        if new_integration:
            new_integrals = self.get_integrals2(codename, *fit_params)
        else:
            new_integrals = self.get_integrals(codename, *fit_params)

        new_heights, new_bins = np.histogram(new_integrals, bins = 50, range = [0,1], normed = True)

        fig = plt.figure(figsize = (16,8))
        sp1 = fig.add_subplot(121)
        sp2 = fig.add_subplot(122)

        sp1.step(bins[:-1] + 0.5*(bins[1]-bins[0]), heights, color = 'k')
        sp1.plot([0,1],[1,1], color = 'r')
        sp1.set_xlim(0,1)
        sp1.set_ylim(0,5)

        sp2.step(new_bins[:-1] + 0.5*(new_bins[1] - new_bins[0]), new_heights, color = 'k')
        sp2.plot([0,1],[1,1], color = 'r')
        sp2.set_xlim(0,1)
        sp2.set_ylim(0,5)

        sp1.set_ylabel('Norm Frequency', family = 'Roboto', weight = 'light', fontsize = 20)

        if new_integration:
            sp2.text(0.02, 0.98, 'Power: %.3f\nBias: %.3f\nShift: %.3f' % fit_params, family = 'Roboto', weight = 'light', fontsize = 20, ha = 'left', va = 'top', transform = sp2.transAxes)
        else:
            sp2.text(0.02, 0.98, 'Stretch: %.3f\nShift: %.3f' % fit_params, family = 'Roboto', weight = 'light', fontsize = 20, ha = 'left', va = 'top', transform = sp2.transAxes)
        fig.text(0.5,0.98, codename.capitalize(), ha = 'center', va = 'center', family = 'Roboto', weight = 'light', fontsize = 24)
        fig.text(0.5, 0.02, 'PIT', family = 'Roboto', weight = 'light', fontsize = 20)

        sp2.set_yticklabels([])

        plt.subplots_adjust(wspace = 0)



    def plot_nz(self, codename, stretch = 1., shift = 0., n_z_bins = np.arange(0, 7, 0.1), chunksize = 2*10**5):

        # NEEDS TO BE UPDATED

        if codename == 'ephor':
            z_mc = self.hsc_cat.ephor_mc
            pdf = self.ephor_pdf
            pdf_bins = self.ephor_bins
        elif codename == 'ephor_ab':
            z_mc = self.hsc_cat.ephor_ab_mc
            pdf = self.ephor_ab_pdf
            pdf_bins = self.ephor_ab_bins
        elif codename == 'demp':
            z_mc = self.hsc_cat.demp_mc
            pdf = self.demp_pdf
            pdf_bins = self.demp_bins
        elif codename == 'frankenz':
            z_mc = self.hsc_cat.frankenz_mc
            pdf = self.frankenz_pdf
            pdf_bins = self.frankenz_bins
        # elif codename == 'mizuki':
        #     z_mc = hsc_cat.mizuki_mc
        # elif codename == 'mlz':
        #     z_mc = hsc_cat.mlz_mc
        # elif codename == 'nnpz':
        #     z_mc = hsc_cat.nnpz_mc

        n_z = np.zeros(len(n_z_bins)-1)


        for x in tqdm(xrange(0, len(pdf), chunksize)):

            chunk_z = z_mc[x:x+chunksize]
            chunk_pdf = pdf[x:x+chunksize]
            chunk_bins = np.array([pdf_bins,]*len(chunk_pdf))

            if stretch != 1:
                chunk_bins = (((chunk_bins.T - chunk_z) * stretch) + chunk_z).T
            if shift != 0:
                chunk_bins = chunk_bins + shift

            chunk_pdf = chunk_pdf[chunk_bins > 0]
            chunk_bins = chunk_bins[chunk_bins > 0]

            norms = np.trapz(chunk_pdf, x = chunk_bins)

            chunk_pdf = chunk_pdf / norms

            #NEEDS FIXING TO ADD PDFS TO N(Z)

            for y in tqdm(xrange(len(n_z_bins) - 1)):

                n_z[y] = n_z[y] + sum(chunk_pdf[(chunk_bins > n_z_bins[y]) & (chunk_bins < n_z_bins[y+1])])

        n_z = n_z * (n_z_bins[1]-n_z_bins[0])**-1.

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.step(n_z_bins[:-1], n_z, color = 'k', linewidth = 2)

        sp.set_title(codename.capitalize(), family = 'Roboto', weight = 'light', fontsize = 26)

        sp.set_xlabel('Redshift (z)', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('N(z)', family = 'Roboto', weight = 'light', fontsize = 24)


    def plot_all_n_z(self, codenames = ['ephor_ab', 'frankenz'], n_z_bins = np.arange(0,7,0.1)):

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        colors = ['#1f77b4', '#ff7f0e']

        for z in tqdm(xrange(len(codenames))):

            n_z_bins, n_z = self.get_n_z(codenames[z], fit = False, n_z_bins = n_z_bins)
            sp.step(n_z_bins[:-1], n_z, color = colors[z], linewidth = 2, linestyle = '--', label = codenames[z].capitalize())

            n_z_bins, n_z = self.get_n_z(codenames[z], fit = True, n_z_bins = n_z_bins)
            sp.step(n_z_bins[:-1], n_z, color = colors[z], linewidth = 2, label = codenames[z].capitalize() + '+PIT')

        sp.legend(loc = 'upper right')

        sp.set_xlabel('Redshift (z)', family = 'Roboto', weight = 'light', fontsize = 24)
        sp.set_ylabel('N(z)', family = 'Roboto', weight = 'light', fontsize = 24)




    def get_n_z(self, codename, fit = False, foldername = './data/n_z_cache/', n_z_bins = np.arange(0,7,0.1), chunksize = 2*10**5):

        if codename == 'ephor_ab':
            z_mc = self.hsc_cat.ephor_ab_mc
            pdf = self.ephor_ab_pdf
            pdf_bins = self.ephor_ab_bins
            fit_params = ephor_ab_PIT_fit
        elif codename == 'frankenz':
            z_mc = self.hsc_cat.frankenz_mc
            pdf = self.frankenz_pdf
            pdf_bins = self.frankenz_bins
            fit_params = frankenz_PIT_fit
        elif codename == 'nnpz':
            z_mc = self.hsc_cat.nnpz_mc
            pdf = self.nnpz_pdf
            pdf_bins = self.nnpz_bins
            fit_params = nnpz_PIT_fit

        if fit:
            fittxt = '_fit'
        else:
            fittxt = ''

        filename = foldername + codename + fittxt + '.dat'

        if not os.path.isdir(foldername):
            os.mkdir(foldername)

        if os.path.isfile(filename):
            # Read in the file
            n_z_bins, n_z = np.loadtxt(filename, unpack = True)
            n_z = n_z[:-1]

        else:

            n_z = np.zeros(len(n_z_bins))
            n_z[-1] = np.nan
            
            if fit:

                stretch, shift = fit_params

                for x in tqdm(xrange(0, len(pdf), chunksize)):

                    chunk_z = z_mc[x:x+chunksize]
                    chunk_pdf = pdf[x:x+chunksize]
                    chunk_bins = np.array([pdf_bins,]*len(chunk_pdf))

                    chunk_bins = (((chunk_bins.T - chunk_z) * stretch) + chunk_z).T
                    chunk_bins = chunk_bins + shift

                    chunk_bins[chunk_bins < 0] = 0

                    chunk_bin_size = chunk_bins[0][-1] - chunk_bins[0][-2]

                    norms = np.trapz(chunk_pdf, x = chunk_bins)

                    chunk_pdf = (chunk_pdf.T / norms).T

                    for y in tqdm(xrange(len(n_z_bins) - 1)):

                        n_z[y] = n_z[y] + sum(chunk_bin_size * chunk_pdf[(chunk_bins >= n_z_bins[y]) & (chunk_bins < n_z_bins[y+1])])

                n_z = n_z * (n_z_bins[1]-n_z_bins[0])**-1.

                np.savetxt(filename, np.vstack((n_z_bins, n_z)).T, fmt = '  %.4e')
                n_z = n_z[:-1]

            else:

                for x in tqdm(xrange(0, len(pdf), chunksize)):

                    chunk_z = z_mc[x:x+chunksize]
                    chunk_pdf = pdf[x:x+chunksize]
                    chunk_bins = np.array([pdf_bins,]*len(chunk_pdf))

                    chunk_bin_size = chunk_bins[0][-1] - chunk_bins[0][-2]

                    norms = np.trapz(chunk_pdf, x = chunk_bins)

                    chunk_pdf = (chunk_pdf.T / norms).T

                    for y in tqdm(xrange(len(n_z_bins) - 1)):

                        n_z[y] = n_z[y] + sum(chunk_bin_size * chunk_pdf[(chunk_bins >= n_z_bins[y]) & (chunk_bins < n_z_bins[y+1])])

                n_z = n_z * (n_z_bins[1]-n_z_bins[0])**-1.

                np.savetxt(filename, np.vstack((n_z_bins, n_z)).T, fmt = '  %.4e')
                n_z = n_z[:-1]

                

        return n_z_bins, n_z