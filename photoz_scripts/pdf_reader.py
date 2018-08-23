from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits
from tqdm import tqdm
import photoz
import os
from scipy.optimize import minimize

hsc_cat = photoz.hsc_cat


class reader:

    # def __init__(self, name = 'ephor'):

    #     self.name = name

    #     inputfp = './data/pdr1_' + name + '_deep_cosmos/'

    #     file_list = sorted(glob(inputfp + '*'))

    #     self.object_id = np.array([])
    #     self.pdf = []
    #     self.bins = np.array([])

    #     idlist = []
    #     pdflist = []
    #     binlist = []

    #     for filename in tqdm(file_list):

    #         hdulist = fits.open(filename)

    #         idlist = idlist + list(hdulist[1].data['ID'])
    #         pdflist.append(hdulist[1].data['PDF'])
    #         binlist.append(hdulist[2].data['BINS'])

    #     if all([all(binlist[0] == rest) for rest in binlist]):
    #         self.bins = binlist[0]

    #     self.pdf = np.vstack(pdflist)

    #     self.object_id = np.array(idlist)


    def __init__(self, name = 'ephor', cachefolder = './data/cache/'):

        self.name = name

        inputfp = './data/pdr1_' + name + '_deep_cosmos/'

        file_list = sorted(glob(inputfp + '*'))

        self.object_id = np.array([])
        self.pdf = []
        self.bins = np.array([])

        idlist = []
        pdflist = []
        binlist = []
        hsc_indices = []

        # This loop goes through each of the files in the PDF directory

        for filename in tqdm(file_list):

            temp_hsc_indices = []
            temp_match_indices = []

            hdulist = fits.open(filename)

            # Set the cache file name for every file in the PDF directory

            cachefile = cachefolder + self.name + '_' + filename.split('/')[-1]
            cachefile = cachefile[:-5] + '.dat'

            # If the cache file doesn't exist, create it

            if not os.path.isfile(cachefile):
                for x in tqdm(xrange(len(hsc_cat.object_id))):

                    # Find which IDs in the PDF file match HSC IDs

                    index = np.where(hsc_cat.object_id[x] == hdulist[1].data['ID'])[0]

                    # If there is an ID match, add it to hsc_indices and match_indices

                    if len(index) != 0:
                        index = index[0]
                        temp_hsc_indices.append(x)
                        temp_match_indices.append(index)

                tabledata = np.array((temp_hsc_indices, temp_match_indices)).T
                if not os.path.isdir('./data/cache/'):
                    os.mkdir('./data/cache/')
                np.savetxt(cachefile, tabledata, header = 'hsc_index pdf_index \n', fmt = '  %8i')

            else:
                temp_hsc_indices, temp_match_indices = np.loadtxt(cachefile, unpack = True, dtype = int)

            hsc_indices = hsc_indices + list(temp_hsc_indices)

            # Extract only the PDFs which have matching IDs with HSC

            pdflist.append(hdulist[1].data['PDF'][temp_match_indices])
            binlist.append(hdulist[2].data['BINS'])

            if all([all(binlist[0] == rest) for rest in binlist]):
                self.bins = binlist[0]

        # Store an attribute which says which objects in the PDF catalog correspond with what objects in the HSC catalog

        self.hsc_indices = np.array(hsc_indices)

        # Turn self.pdf in to an array composed of many PDFs along axis 0

        self.pdf = np.vstack(pdflist)



    def get_integrals(self, stretch = 1., shift = 0.):

        if self.name == 'ephor':
            z_mc = hsc_cat.ephor_mc
        elif self.name == 'ephor_ab':
            z_mc = hsc_cat.ephor_ab_mc
        elif self.name == 'demp':
            z_mc = hsc_cat.demp_mc
        elif self.name == 'frankenz':
            z_mc = hsc_cat.frankenz_mc
        elif self.name == 'mizuki':
            z_mc = hsc_cat.mizuki_mc
        elif self.name == 'mlz':
            z_mc = hsc_cat.mlz_mc
        elif self.name == 'nnpz':
            z_mc = hsc_cat.nnpz_mc

        integrals = []

        # Calculate the cumulative integral from the lower bound to the estimated redshift

        for x in xrange(len(self.pdf)):

            mod_pdf = self.pdf[x]
            mod_bins = self.bins
            this_montecarlo = z_mc[self.hsc_indices[x]]

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



    def plot_PIT(self, shift = 0., stretch = 1.):

        integrals = self.get_integrals(stretch, shift)

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        sp.hist(integrals, bins = 50, normed = True, range = [0,1], histtype = 'step', color = 'k')
        sp.plot([0,1],[1,1], color = 'r')
        sp.set_xlim(0,1)
        sp.set_ylim(0,5)



    def plot_PIT_fit(self):

        fig = plt.figure(figsize = (16,8))
        sp1 = fig.add_subplot(121)
        sp2 = fig.add_subplot(122)

        integrals = self.get_integrals()

        heights, bins = np.histogram(integrals, bins = 50, range = [0,1], normed = True)

        sp1.step(bins[:-1] + 0.5*(bins[1]-bins[0]), heights, color = 'k')
        sp1.plot([0,1],[1,1], color = 'r')
        sp1.set_xlim(0,1)
        sp1.set_ylim(0,5)


        def calchist((stretch, shift)):

            mod_integrals = self.get_integrals(stretch, shift)
            mod_heights, mod_bins = np.histogram(mod_integrals, bins = 50, range = [0,1], normed = True)

            return sum((mod_heights-1.)**2.)

        result = minimize(calchist, [1,0], method = 'Nelder-Mead')

        fit_stretch, fit_shift = result.x

        print result

        new_integrals = get_integrals(fit_stretch, fit_shift)

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


