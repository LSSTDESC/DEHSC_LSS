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

    def __init__(self, name = 'ephor'):

        self.name = name

        inputfp = './data/pdr1_' + name + '_deep_cosmos/'

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



def get_matches(pdfreader, cachefile = './data/cache/PIT_cache.dat'):

    if os.path.isfile(cachefile):
        hsc_indices, match_indices = np.loadtxt(cachefile, unpack = True)

    else:

        hsc_indices = [] # The HSC objects with matches in the pdf indices
        match_indices = [] # The matching pdf indices

        for x in tqdm(xrange(len(hsc_cat.object_id))):

            index = np.where(hsc_cat.object_id[x] == pdfreader.object_id)[0]

            if len(index) != 0:
                index = index[0]
                hsc_indices.append(x)

                match_indices.append(index)

        # mod_bins = [pdfreader.bins,] * len(hsc_indices)
        # mod_pdfs = []

        # for x in tqdm(xrange(len(hsc_indices))):

        #     mod_bins[x] = mod_bins[x][mod_bins[x] < hsc_cat.ephor_median[hsc_indices[x]]]
        #     mod_pdfs.append(pdfreader.pdf[match_indices[x]][:len(mod_bins[x])])

        # mod_pdfs = np.array(mod_pdfs)
        # mod_bins = np.array(mod_bins)
        match_indices = np.array(match_indices)

        tabledata = np.array((hsc_indices, match_indices)).T
        if not os.path.isdir('./data/cache/'):
            os.mkdir('./data/cache/')
        np.savetxt(cachefile, tabledata, header = 'hsc_index pdf_index integrals \n', fmt = '  %8i  %8i')

    return hsc_indices, match_indices



def get_integrals(hsc_indices, match_indices, pdfreader, stretch = 1., shift = 1.):

    if pdfreader.name == 'ephor':
        z_mc = hsc_cat.ephor_mc
    elif pdfreader.name == 'ephor_ab':
        z_mc = hsc_cat.ephor_ab_mc
    elif pdfreader.name == 'demp':
        z_mc = hsc_cat.demp_mc
    elif pdfreader.name == 'frankenz':
        z_mc = hsc_cat.frankenz_mc
    elif pdfreader.name == 'mizuki':
        z_mc = hsc_cat.mizuki_mc
    elif pdfreader.name == 'mlz':
        z_mc = hsc_cat.mlz_mc
    elif pdfreader.name == 'nnpz':
        z_mc = hsc_cat.nnpz_mc

    integrals = []

    # Calculate the cumulative integral from the lower bound to the estimated redshift

    for x in xrange(len(match_indices)):

        mod_pdf = pdfreader.pdf[match_indices[x]]
        mod_bins = pdfreader.bins

        if stretch != 1.:
            mod_bins = ((mod_bins - z_mc[hsc_indices[x]]) * stretch) + z_mc[hsc_indices[x]]
            mod_pdf = mod_pdf[mod_bins > 0]
            mod_bins = mod_bins[mod_bins > 0]

        if shift != 1.:
            mod_bins = mod_bins + shift

        part = np.trapz(mod_pdf[pdfreader.bins < z_mc[hsc_indices[x]]], x = mod_bins[pdfreader.bins < z_mc[hsc_indices[x]]])
        total = np.trapz(mod_pdf, x = mod_bins)
        integrals.append(part/total)

    integrals = np.array(integrals)

    return integrals



def plot_PIT(pdfreader):

    hsc_indices, match_indices = get_matches(pdfreader)
    integrals = get_integrals(hsc_indices, match_indices, pdfreader)

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    sp.hist(integrals, bins = 50, normed = True, range = [0,1], histtype = 'step', color = 'k')
    sp.plot([0,1],[1,1], color = 'r')
    sp.set_xlim(0,1)
    sp.set_ylim(0,5)



def plot_PIT_fit(pdfreader):

    fig = plt.figure(figsize = (16,8))
    sp1 = fig.add_subplot(121)
    sp2 = fig.add_subplot(122)

    hsc_indices, match_indices = get_matches(pdfreader)
    integrals = get_integrals(hsc_indices, match_indices, pdfreader)

    heights, bins = np.histogram(integrals, bins = 50, range = [0,1], normed = True)

    sp1.step(bins[-1] + 0.5*(bins[1]-bins[0]), heights, color = 'k')
    sp1.plot([0,1],[1,1], color = 'r')
    sp1.set_xlim(0,1)
    sp1.set_ylim(0,5)


    def calchist(stretch, shift):

        mod_integrals = get_integrals(hsc_indices, match_indices, pdfreader, stretch, shift)
        mod_heights, mod_bins = np.histogram(mod_integrals, bins = 50, range = [0,1], normed = True)

        return sum((mod_heights-1)**2)

    fit_stretch, fit_shift = minimize(calchist, [1,1]).x

    new_integrals = get_integrals(hsc_indices, match_indices, pdfreader, fit_stretch, fit_shift)

    new_heights, new_bins = np.histogram(new_integrals, bins = 50, range = [0,1], normed = True)

    sp2.step(bins[-1] + 0.5*(bins[1] - bins[0]), heights, color = 'k')
    sp2.plot([0,1],[1,1], color = 'r')
    sp2.set_xlim(0,1)
    sp2.set_ylim(0,5)



