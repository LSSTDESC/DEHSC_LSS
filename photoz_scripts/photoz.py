from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import fits
from os.path import isfile
from tqdm import tqdm
from scipy.optimize import curve_fit

font = {'family':'Roboto', 'weight':'light'}





def build_cosmos_cat(inputfp = './data/COSMOS2015_Laigle+_v1.1.fits', outputfp = './data/photoz_COSMOS.dat', hsc_paper_cuts = True):

    """
        This function reads in the full COSMOS catalog from inputfp and rewrites it in a more easily usable form at outputfp.
        Additionally, it excludes stars (z=0), X-ray detected sources (z=9) and objects with poor fits (z<0).
    """


    hdulist = fits.open(inputfp)
    data = hdulist[1].data

    ID = data['NUMBER']
    ra = data['ALPHA_J2000']
    dec = data['DELTA_J2000']
    rmag = data['r_MAG_AUTO']
    photoz = data['PHOTOZ']

    if hsc_paper_cuts == True:
        # See https://arxiv.org/pdf/1704.05988.pdf 3rd page on "COSMOS data" part 2, "Photo-z"
        lim_indices = np.where((0.01 < photoz) & (9 > photoz) & (data['TYPE'] == 0) & (data['ZP_2'] < 0) & (data['MASS_BEST'] > 7.5) 
            & (np.maximum(data['ZPDF_H68'] - data['ZPDF'], data['ZPDF'] - data['ZPDF_L68']) < 0.05*(1+photoz))
            & (data['CHI2_BEST'] < data['CHIS']) & (data['CHI2_BEST']/data['NBFILT'] < 5.))

        ID = ID[lim_indices]
        ra = ra[lim_indices]
        dec = dec[lim_indices]
        rmag = rmag[lim_indices]
        photoz = photoz[lim_indices]

    tabledata = np.array((ID, ra, dec, rmag, photoz)).T
    np.savetxt(outputfp, tabledata, header = 'id ra dec rmag photoz\n', fmt = '  %8i  %8f  %8f  %5f  %5f')



class cosmos_reader:

    """
        Generates the file at inputfp if it does not exist, or otherwise reads it in.  Object attributes are named according to their columns.
    """

    def __init__(self, inputfp = './data/photoz_COSMOS.dat'):

        if not isfile(inputfp):
            build_cosmos_cat()

        self.id, self.ra, self.dec, self.mag_r, self.z_phot = pd.read_csv(inputfp, header = None, comment = '#', delimiter = '\s+').values.T
        self.valid_index = np.where(self.z_phot > 0)[0]


class hsc_reader:

    """
        Reads in the data from an HSC csv file generated from a SQL query.  Values can be added to this object, but be careful about taking them
        out because lots of the other methods require the values already put here.
    """

    def __init__(self, inputfp = './data/11473.csv'):

        (self.object_id, self.specz_id, self.ra, self.dec, 
            self.mag_g, self.mag_r, self.mag_i, self.mag_z, self.mag_y, 
            self.z_spec, self.ephor_mc, self.ephor_ab_mc, self.demp_mc,
            self.frankenz_mc, self.mizuki_mc, self.mlz_mc, self.nnpz_mc,
            self.flag_zcosmos, self.flag_cosmos_fmos) = pd.read_csv(inputfp).values.T


hsc_cat = hsc_reader()
cosmos_cat = cosmos_reader()


def match_cat():

    """
        Matches indices between the HSC catalog and the COSMOS photometric catalog.  The first returned item is the indices in the cosmos_reader object
        which will match them to each object in the HSC catalog.  The second returned item is the distance between each pair of matched objects in 
        arcseconds.
    """

    #Force to match only to galaxies

    gals = np.where((cosmos_cat.z_phot > 0.001) & (cosmos_cat.z_phot < 9.9))[0]
    cosmos_skycoord = SkyCoord(ra = cosmos_cat.ra[gals] * u.deg, dec = cosmos_cat.dec[gals] * u.deg)

    hsc_skycoord = SkyCoord(ra = hsc_cat.ra * u.deg, dec = hsc_cat.dec * u.deg)
    
    cosmos_index, dist_2d, dist_3d = hsc_skycoord.match_to_catalog_sky(cosmos_skycoord)

    return gals[cosmos_index], dist_2d.to('arcsec').value




def dist_hist():

    """
        Plot a histogram of all of the (logged) distances between matched galaxies.
    """

    cosmos_index, dist_2d = match_cat()

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    sp.hist(np.log10(dist_2d), bins = 30, histtype = 'step', color = 'k', linewidth = 2)

    sp.set_xlabel('$log_{10}$[Separation/Arcsec]', fontdict = font, fontsize = 24)
    sp.set_ylabel('Frequency', fontdict = font, fontsize = 24)



def sep_mag_diff():

    """
        Plot the 2d separation between galaxy matches and the difference in their magnitudes in r-band.
    """
    
    cosmos_index, dist_2d = match_cat()

    fig = plt.figure(figsize = (8,8))
    sp = fig.add_subplot(111)

    cosmos_r = cosmos_cat.mag_r[cosmos_index]
    hsc_r = hsc_cat.mag_r

    # Make sure r-band magnitude is a valid value

    dist_2d = dist_2d[cosmos_r < 50]
    hsc_r = hsc_r[cosmos_r < 50]
    cosmos_r = cosmos_r[cosmos_r < 50]

    magdiff = hsc_r - cosmos_r

    sp.scatter(np.log10(dist_2d), magdiff, edgecolors = 'None', facecolors = 'k', marker = '.')

    sp.set_ylabel('$r_{HSC} - r_{COSMOS}$', fontdict = font, fontsize = 24)
    sp.set_xlabel('$log_{10}$[Separation/Arcsec]', fontdict = font, fontsize = 24)






def redshift_compare_grid(rmaglims = [24,25], seplims = [1,2,4]):

    rmaglims = [-np.inf] + rmaglims + [np.inf]
    seplims = [-np.inf] + seplims + [np.inf]

    fig, subplots = plt.subplots(len(rmaglims)-1, len(seplims)-1, figsize = (4*(len(rmaglims)-1), 4*(len(seplims))-1))

    cosmos_index, dist_2d = match_cat()

    for x in tqdm(xrange(len(rmaglims)-1)):

        for y in tqdm(xrange(len(seplims)-1)):

            lim_indices = np.where((cosmos_cat.mag_r[cosmos_index] > rmaglims[x]) & (cosmos_cat.mag_r[cosmos_index] < rmaglims[x+1]) &
             (dist_2d > seplims[y]) & (dist_2d < seplims[y+1]))

            hsc_redshift = hsc_cat.z_spec[lim_indices]
            cosmos_redshift = cosmos_cat.z_phot[cosmos_index][lim_indices]

            subplots[x][y].scatter(hsc_redshift, cosmos_redshift, edgecolors = 'None', facecolors = 'k', marker = '.')
            subplots[x][y].set_xlim(-0.5, 10.5)
            subplots[x][y].set_ylim(-0.5, 10.5)

            f_in = float(len(np.where(abs(hsc_redshift - cosmos_redshift) < 0.1*(1+hsc_redshift))[0]))/float(len(hsc_redshift))

            subplots[x][y].text(0.98, 0.98, r'$f_{\Delta z < 0.1} =' + ' %.2f$' % f_in, fontdict = font, fontsize = 15, ha = 'right', va = 'top', transform = subplots[x][y].transAxes)

            if x != len(rmaglims) - 2:
                subplots[x][y].set_xticklabels([])
            if y != 0:
                subplots[x][y].set_yticklabels([])

    plt.subplots_adjust(hspace = 0, wspace = 0)

    subplots[0][0].set_ylabel('r < 24', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[1][0].set_ylabel('24 < r < 25', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[2][0].set_ylabel('r > 25', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[0][0].set_title('Sep < 1"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][1].set_title('1" < Sep < 2"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][2].set_title('2" < Sep < 4"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][3].set_title('Sep > 4"', family = 'Roboto', weight = 'light', fontsize = 20)

    fig.text(0.5, 0.03, 'HSC specz', family = 'Roboto', weight = 'light', fontsize = 24)
    fig.text(0.03, 0.5, 'COSMOS photoz', fontsize = 24, fontdict = font, ha = 'center', va = 'center', rotation = 'vertical')



def delta_redshift_compare_grid(rmaglims = [24,25], seplims = [1,2,4], fit = True):

    rmaglims = [-np.inf] + rmaglims + [np.inf]
    seplims = [-np.inf] + seplims + [np.inf]

    fig, subplots = plt.subplots(len(rmaglims)-1, len(seplims)-1, figsize = (4*(len(rmaglims)-1), 4*(len(seplims))-1))

    cosmos_index, dist_2d = match_cat()

    for x in tqdm(xrange(len(rmaglims)-1)):

        for y in tqdm(xrange(len(seplims)-1)):

            lim_indices = np.where((cosmos_cat.mag_r[cosmos_index] > rmaglims[x]) & (cosmos_cat.mag_r[cosmos_index] < rmaglims[x+1]) &
             (dist_2d > seplims[y]) & (dist_2d < seplims[y+1]))

            hsc_redshift = hsc_cat.z_spec[lim_indices]
            cosmos_redshift = cosmos_cat.z_phot[cosmos_index][lim_indices]

            bias = (hsc_redshift - cosmos_redshift)/(1+hsc_redshift)

            heights, bins = np.histogram(bias, bins = 30, range = [-.2,.2], normed = True)
            binctr = bins[:-1]+ 0.5*(bins[1]-bins[0])

            subplots[x][y].step(binctr, heights, color = 'b', linewidth = 2)
            subplots[x][y].set_xlim(-0.21, 0.21)
            subplots[x][y].set_ylim(0.01, 35)

            # f_in = float(len(np.where(abs(hsc_redshift - cosmos_redshift) < 0.1*(1+hsc_redshift))[0]))/float(len(hsc_redshift))
            # sigma = np.std(bias[(bias > -0.6) & (bias < 0.6)])
            sigma = np.std(bias)

            if fit:
                def gaussfit(x,a,b):
                    return a * (b*np.sqrt(2*np.pi))**-1 * np.exp(-0.5*x**2/b**2)

                params, cov = curve_fit(gaussfit, binctr, heights)
                xdata = np.linspace(-.2,.2,1000)
                subplots[x][y].plot(xdata, gaussfit(xdata, *params), color = 'r', label = r'$\sigma = ' + '%.3f$' % params[1])
                subplots[x][y].legend(loc = 'upper left', fontsize = 10)

            else:            
                subplots[x][y].text(0.98, 0.98, r'$\sigma = ' + '%.2f$' % sigma, fontdict = font, fontsize = 15, ha = 'right', va = 'top', transform = subplots[x][y].transAxes)

            if x != len(rmaglims) - 2:
                subplots[x][y].set_xticklabels([])
            if y != 0:
                subplots[x][y].set_yticklabels([])

    plt.subplots_adjust(hspace = 0, wspace = 0)

    subplots[0][0].set_ylabel('r < 24', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[1][0].set_ylabel('24 < r < 25', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[2][0].set_ylabel('r > 25', family = 'Roboto', weight = 'light', fontsize = 20, x = 0.02)
    subplots[0][0].set_title('Sep < 1"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][1].set_title('1" < Sep < 2"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][2].set_title('2" < Sep < 4"', family = 'Roboto', weight = 'light', fontsize = 20)
    subplots[0][3].set_title('Sep > 4"', family = 'Roboto', weight = 'light', fontsize = 20)

    fig.text(0.5, 0.03, 'Redshift Error', family = 'Roboto', weight = 'light', fontsize = 24)
    fig.text(0.03, 0.5, 'Normalized Frequency', fontsize = 24, fontdict = font, ha = 'center', va = 'center', rotation = 'vertical')
