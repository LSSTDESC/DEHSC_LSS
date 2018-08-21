from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from glob import glob
from astropy.io import fits
from tqdm import tqdm
import photoz
import os

hsc_cat = photoz.hsc_cat


class reader:

    def __init__(self, inputfp = './data/pdr1_ephor_deep_cosmos/'):

        file_list = sorted(glob(inputfp + '*'))

        self.object_id = np.array([])
        self.pdf = []
        self.bins = np.array([])

        idlist = []
        pdflist = []
        binlist = []

        for filename in tqdm(file_list):

            hdulist = fits.open(filename)

            idlist = idlist + list(hdulist[1].data['ID'])
            pdflist.append(hdulist[1].data['PDF'])
            binlist.append(hdulist[2].data['BINS'])

        if all([all(binlist[0] == rest) for rest in binlist]):
            self.bins = binlist[0]

        self.pdf = np.concatenate(tuple(pdflist))

        self.object_id = np.array(idlist)



def get_PIT(pdfreader, cachefile = './data/cache/PIT_cache.dat'):

    if os.path.isfile(cachefile):
        hsc_indices, match_indices, integrals = np.loadtxt(cachefile, unpack = True)

    else:

        hsc_indices = [] # The HSC objects with matches in the pdf indices
        match_indices = [] # The matching pdf indices

        for x in tqdm(xrange(len(hsc_cat.ephor_mc))):

            index = np.where(hsc_cat.object_id[x] == pdfreader.object_id)[0]

            if len(index) != 0:
                index = index[0]
                hsc_indices.append(x)

                match_indices.append(index)

        mod_bins = [pdfreader.bins,] * len(hsc_indices)
        mod_pdfs = []

        for x in tqdm(xrange(len(hsc_indices))):

            mod_bins[x] = mod_bins[x][mod_bins[x] < hsc_cat.ephor_median[hsc_indices[x]]]
            mod_pdfs.append(pdfreader.pdf[index][:len(mod_bins[x])])

        mod_pdfs = np.array(mod_pdfs)
        match_indices = np.array(match_indices)
        mod_bins = np.array(mod_bins)

        integrals = []

        # Calculate the cumulative integral from the lower bound to the estimated redshift

        # integrals = np.trapz(mod_pdfs, x = mod_bins)/np.trapz(pdfreader.pdf[match_indices], x = pdfreader.bins[match_indices])

        for x in tqdm(xrange(len(mod_bins))):
            integrals.append(np.trapz(mod_pdfs[x], x = mod_bins[x])/np.trapz(pdfreader.pdf[match_indices[x]], x = pdfreader.bins))

        integrals = np.array(integrals)

        tabledata = np.array((hsc_indices, match_indices, integrals)).T
        os.mkdir('./data/cache/')
        np.savetxt(cachefile, tabledata, header = 'hsc_index pdf_index integrals \n', fmt = '  %8i  %8i  %8f')

    return hsc_indices, match_indices, integrals


def plot_PIT(pdfreader):

    hsc_indices, match_indices, integrals = get_PIT(pdfreader)

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    sp.hist(integrals, bins = 50, normed = True, histtype = 'step', color = 'k')
    sp.plot([0,1],[1,1], color = 'r')
