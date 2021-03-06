# Make sure we can find ceci and NaMaster

import sys
sys.path.append('/global/homes/a/abrouss/ceci')
sys.path.append('/global/homes/a/abrouss/HSCLink')
sys.path.append('/global/u1/a/abrouss/NaMaster')
# sys.path.append('/global/homes/a/abrouss/HyperSuprimeStructure-HSC-LSS/hsc_lss/flatmaps.py')
# sys.path.append('/global/u1/a/abrouss/.local/cori/3.6-anaconda-4.4/lib/python3.6/site-packages')

from glob import glob
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from numpy.random import poisson
from hsc_lss.flatmaps import read_flat_map, FlatMapInfo
from hsc_lss.map_utils import createCountsMap
from hsc_lss.tracer import Tracer
import pymaster as nmt
from tqdm import tqdm

import matplotlib as mpl
from matplotlib import pyplot as plt


catalogs = glob('/global/cscratch1/sd/damonge/HSC_ceci/WIDE_*_sirius_out/')

ell_bins = [100.0,200.0,300.0,400.0,600.0,800.0,1000.0,1400.0,1800.0,2200.0,3000.0,3800.0,4600.0,6200.0,7800.0,9400.0,12600.0,15800.0]
pz_bins = [0.15,0.50,0.75,1.00,1.50]


class catalog:
    def __init__(self, catfolder):
        """
        Holds all of the information relevant to a given HSC field
        """
        self.name = catfolder.split('_')[-3]
        self.filepath = catfolder

        self.cat = fits.open(self.filepath + 'clean_catalog.fits')[1].data

        self.masked_fraction_fp = self.filepath + 'masked_fraction.fits'
        self.masked_fraction = fits.open(self.masked_fraction_fp)[0].data.ravel()
        self.mask_binary = np.ones_like(self.masked_fraction)
        self.mask_binary[self.masked_fraction<0.5] = 0.

        self.weight = self.masked_fraction * self.mask_binary
        self.goodpix = np.where(self.mask_binary>0.1)[0]

        self.fsk,_ = read_flat_map(self.masked_fraction_fp)

        self.nmaps, self.nmaps_s1, self.nmaps_s2 = self.get_nmaps_split()

        self.fields = [field(self, thisbin) for thisbin in self.nmaps]
        self.fields_s1 = [field(self, thisbin) for thisbin in self.nmaps_s1]
        self.fields_s2 = [field(self, thisbin) for thisbin in self.nmaps_s2]

        self.ell_bins = nmt.NmtBinFlat(ell_bins[:-1], ell_bins[1:])
        self.ell_bins_uncpld = self.ell_bins.get_effective_ells()

        self.wsp = nmt.NmtWorkspaceFlat()
        self.wsp.compute_coupling_matrix(self.fields[0].field, self.fields[0].field, self.ell_bins)



    # Modified from catmapper.py
    def get_nmaps_split(self, pz_code = 'ephor_ab', pz_mark = 'best'):
        """
        Calculates the number map for a given catalog and splits it into two halves.

        Returns the full number map, 1st split, and 2nd split
        """

        maps = []
        maps_split1 = []
        maps_split2 = []

        if pz_code == 'ephor_ab':
            pz_code_col = 'eab'

        column_mark='pz_'+pz_mark+'_'+pz_code_col

        for zi,zf in zip(pz_bins[:-1], pz_bins[1:]) :
            msk_bin=(self.cat[column_mark]<=zf) & (self.cat[column_mark]>zi)
            subcat=self.cat[msk_bin]
            index_split1 = np.random.choice(np.arange(len(subcat['ra'])), int(len(subcat['ra'])/2), replace = False)
            index_split2 = np.delete(np.arange(len(subcat['ra'])), index_split1)
            nmap=createCountsMap(subcat['ra'],subcat['dec'],self.fsk)
            nmap_split1 = createCountsMap(subcat['ra'][index_split1], subcat['dec'][index_split1], self.fsk)
            nmap_split2 = createCountsMap(subcat['ra'][index_split2], subcat['dec'][index_split2], self.fsk)
            maps.append(nmap)
            maps_split1.append(nmap_split1)
            maps_split2.append(nmap_split2)
        return np.array(maps), np.array(maps_split1), np.array(maps_split2)



    def get_power_spectra(self, field1, field2):

        """
        Calculates the power spectrum for a given pair of number maps
        """

        c_l = []

        for zbin1, zbin2 in zip(field1, field2):

            coupled_c_l = nmt.compute_coupled_cell_flat(zbin1.field, zbin2.field, self.ell_bins)
            c_l.append(self.wsp.decouple_cell(coupled_c_l)[0])

        return np.array(c_l)




    def get_power_spectra_all(self):

        """
        Calculates the power spectra for various combinations of the original and split number maps
        """

        pair_list = [(self.fields, self.fields), (self.fields_s1, self.fields_s1), 
                (self.fields_s2, self.fields_s2), (self.fields_s1, self.fields_s2)]
        
        cl_list = []

        for thisfield1, thisfield2 in pair_list:

            cl_list.append(self.get_power_spectra(thisfield1, thisfield2))

        return cl_list



    # Modified from powerspecter.py
    def calc_shotnoise(self, fieldlist):
        """
        Calculates the theoretical shot noise for a list of redshift bins
        """

        nbins = len(pz_bins)-1
        nell = len(self.ell_bins_uncpld)

        nls_all=np.zeros([nbins, nell])
        for i in range(nbins) :
            corrfac=np.sum(self.weight)/(self.fsk.nx*self.fsk.ny)
            nl=np.ones(nell)*corrfac/fieldlist[i].ndens_perad
            nls_all[i]=self.wsp.decouple_cell([nl])[0]
        return nls_all



    def plot_power_spectra(self, yaxis = 'log'):

        """
        Plots the power spectra for <N,N>, <s1,s1>, <s2,s2>, <s1,s2> where N is the original number map and s1, s2 are the two split halves
        """

        cl_list = self.get_power_spectra_all()

        fig, subplots = plt.subplots(2,2, figsize = (8,8))

        subplots = subplots.flatten()
        names = ['0,0', 's1,s1', 's2,s2', 's1,s2']

        for x, (sp, fieldcl, name) in enumerate(zip(subplots, cl_list, names)):
            for y, zbin_cl in enumerate(fieldcl):
                sp.plot(self.ell_bins_uncpld, zbin_cl)
                sp.text(0.02, 0.02, name, transform = sp.transAxes, fontsize = 18, ha = 'left', va = 'bottom')

            sp.set_xscale('log')
            sp.set_yscale(yaxis)
            # sp.set_xlabel(r'$\ell$')
            # sp.set_ylabel(r'C_$\ell$')


        fig.text(0.5, 0.05, r'$\ell$', fontsize = 24)
        fig.text(0.06, 0.5, r'C$_\ell$', fontsize = 24, ha = 'center', va = 'center', rotation = 'vertical')
        fig.text(0.5, 0.9, self.name, fontsize = 24, ha = 'center', va = 'bottom')

        return fig, subplots



    def plot_shotnoise_test(self):

        """
        Plots the shot noise estimate from 0.25*(<s1,s1> + <s2,s2> - 2<s1,s2>) and the theoretically calculated shot noise for a catalog
        """

        cl_list = self.get_power_spectra_all()

        fig = plt.figure(figsize = (8,8))
        sp = fig.add_subplot(111)

        s1s1 = cl_list[1]
        s2s2 = cl_list[2]
        s1s2 = cl_list[3]

        for x, (zbin11, zbin22, zbin12, noise) in enumerate(zip(s1s1, s2s2, s1s2, self.calc_shotnoise(self.fields))):

            sp.plot(self.ell_bins_uncpld, noise, linewidth = 1.5, color = 'C' + str(x), zorder = 0, alpha = 0.4)
            sp.plot(self.ell_bins_uncpld, 0.25 * (zbin11 + zbin22 - 2.*zbin12), linewidth = 2, color = 'C' + str(x), zorder = 1)

        sp.set_xscale('log')
        sp.set_yscale('log')

        sp.plot([],[],linewidth = 1.5, color = 'C0', alpha = 0.4, label = 'Analytic')
        sp.plot([],[],linewidth = 2, color = 'C0', label = 'Split Map')

        sp.legend(loc = 'upper right', fontsize = 18)

        sp.set_xlabel(r'$\ell$', fontsize = 20)
        sp.set_ylabel('Estimated Shot Noise', fontsize = 20)
        sp.set_title(self.name, fontsize = 24)

        # fig.text(0.5, 0.03, r'$\ell$', fontsize = 20)
        # fig.text(0.06, 0.5, r'<s1,s1> + <s2, s2> - 2 <s1, s2>', fontsize = 20, ha = 'center', va = 'center', rotation = 'vertical')
        # fig.text(0.02, 0.5, r'Estimated Shot Noise', fontsize = 20, ha = 'center', va = 'center', rotation = 'vertical')
        # fig.text(0.5, 0.9, catname, fontsize = 24, ha = 'center', va = 'bottom')

        return fig, sp




class field:
    def __init__(self, catalog, nmap_zbin):
        
        """
        Holds all of the relevant information for a given redshift bin number map
        """

        Ngal = np.sum(nmap_zbin * catalog.mask_binary)
        self.ndens = Ngal/np.sum(catalog.weight)
        self.ndens_perad=self.ndens/(np.radians(catalog.fsk.dx)*np.radians(catalog.fsk.dy))
        self.delta = np.zeros_like(catalog.weight)
        self.delta[catalog.goodpix] = nmap_zbin[catalog.goodpix]/(self.ndens*catalog.masked_fraction[catalog.goodpix])-1

        self.field = nmt.NmtFieldFlat(np.radians(catalog.fsk.lx),np.radians(catalog.fsk.ly),
                                        catalog.weight.reshape([catalog.fsk.ny,catalog.fsk.nx]),
                                        [self.delta.reshape([catalog.fsk.ny,catalog.fsk.nx])])



def plot_shotnoise_avg_test():

    """
    Plots the average shot noise estimates across all of the HSC fields
    """

    ells = None
    theoretical_noise = [[] for x in range(len(pz_bins)-1)]
    shot_estimate = [[] for x in range(len(pz_bins)-1)]

    for catnum, catalog_fp in enumerate(tqdm(catalogs)):

        thiscat = catalog(catalog_fp)

        if not hasattr(ells, '__iter__'):
            ells = thiscat.ell_bins_uncpld

        cl_list = thiscat.get_power_spectra_all()

        s1s1 = cl_list[1]
        s2s2 = cl_list[2]
        s1s2 = cl_list[3]

        for zbin11, zbin22, zbin12, noise, this_theor_noise, this_shot_est in zip(s1s1, s2s2, s1s2, thiscat.calc_shotnoise(thiscat.fields), theoretical_noise, shot_estimate):

            this_theor_noise.append(noise)
            this_shot_est.append(0.25 * (zbin11 + zbin22 - 2.*zbin12))

    theoretical_noise = np.array(theoretical_noise)
    shot_estimate = np.array(shot_estimate)

    avg_theoretical_noise = np.average(theoretical_noise, axis = 1)
    avg_shot_estimate = np.average(shot_estimate, axis = 1)

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    for znum, (this_atn, this_ase) in enumerate(zip(avg_theoretical_noise, avg_shot_estimate)):

        sp.plot(ells, avg_theoretical_noise[znum], linewidth = 1.5, color = 'C' + str(znum), zorder = 0, alpha = 0.4)
        sp.plot(ells, avg_shot_estimate[znum], linewidth = 2, color = 'C' + str(znum), zorder = 1)

    sp.set_xscale('log')
    sp.set_yscale('log')

    sp.plot([],[],linewidth = 1.5, color = 'C0', alpha = 0.4, label = 'Analytic')
    sp.plot([],[],linewidth = 2, color = 'C0', label = 'Split Map')

    sp.legend(loc = 'upper right', fontsize = 18)

    sp.set_xlabel(r'$\ell$', fontsize = 20)
    sp.set_ylabel('Estimated Shot Noise', fontsize = 20)
    
    return fig, sp