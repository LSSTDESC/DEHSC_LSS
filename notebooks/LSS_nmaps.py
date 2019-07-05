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

catalogs = glob('/global/cscratch1/sd/damonge/HSC_ceci/WIDE_*_sirius_out/')

ell_bins = [100.0,200.0,300.0,400.0,600.0,800.0,1000.0,1400.0,1800.0,2200.0,3000.0,3800.0,4600.0,6200.0,7800.0,9400.0,12600.0,15800.0]
z_bins = [0.15,0.50,0.75,1.00,1.50]


# modified from cat_mapper
def get_nmaps(catfolder, pz_code = 'ephor_ab', pz_mark = 'best') :
    """
    Get number counts map from catalog
    """
    maps=[]
    
    if pz_code == 'ephor_ab':
        pz_code_col = 'eab'
    
    cat=fits.open(catfolder + '/clean_catalog.fits')[1].data
    column_mark='pz_'+pz_mark+'_'+pz_code_col

    for zi,zf in zip(pz_i, pz_f) :
        msk_bin=(cat[column_mark]<=zf) & (cat[column_mark]>zi)
        subcat=cat[msk_bin]
        masked_fraction = catfolder + '/masked_fraction.fits'
        fsk,_=read_flat_map(masked_fraction)
        nmap=createCountsMap(subcat['ra'],subcat['dec'],fsk)
        maps.append(nmap)
    return np.array(maps)



def get_nmaps_split(catfolder, pz_bins = [0.15,0.50,0.75,1.00,1.50]):

    pz_i = pz_bins[:-1]
    pz_f = pz_bins[1:]

    if catfolder[-1] != '/':
        catfolder = catfolder + '/'

    thismap = get_nmaps(catfolder)
    masked_fraction = catfolder + 'masked_fraction.fits'
    fsk,_=read_flat_map(masked_fraction)
    geometry = (fsk.ny, fsk.nx)

    newmaps = []
    newmaps_split1 = []
    newmaps_split2 = []

    cutcounter = 0
    bigcutcounter = 0

    for binmap in thismap:
        binmap_split1 = np.zeros(len(binmap))
        
        binmap_split1 = poisson(binmap/2.)
        binmap_split2 = binmap - binmap_split1
        cutcounter += len(np.where((binmap_split1 - binmap) > 0)[0])
        bigcutcounter += len(np.where((binmap_split1 - binmap) > 2)[0])
        
        # Take care of "negative" galaxy count pixels
        # Get rid of the negative galaxy counts while keeping binmap_split1 + binmap_split2 = binmap
        
        negative = np.where(binmap_split2 < 0)
        binmap_split1[negative] = binmap_split1[negative] + binmap_split2[negative]
        binmap_split2[negative] = 0

        newmaps.append(binmap)
        newmaps_split1.append(binmap_split1)
        newmaps_split2.append(binmap_split2)

    return newmaps, newmaps_split1, newmaps_split2



def get_nmtmap(nmap_zbin):

    masked_fraction_file = catfolder + 'masked_fraction.fits'
    masked_fraction = fits.open(masked_fraction_file)[0].data
    fsk,_=read_flat_map(masked_fraction_file)

    mask_binary = np.ones_like(masked_fraction)

    weight = masked_fraction * mask_binary
    goodpix = np.where(mask_binary>0.1)[0]
    Ngal = np.sum(nmap_zbin * mask_binary)
    ndens = Ngal/np.sum(weight)
    delta = np.zeros_like(weight)
    delta[goodpix] = nmap_zbin[goodpix]/(ndens*masked_fraction[goodpix])-1

    field=nmt.NmtFieldFlat(np.radians(fsk.lx),np.radians(fsk.ly),
                                    weight.reshape([fsk.ny,fsk.nx]),
                                    [delta.reshape([fsk.ny,fsk.nx])],templates=conts)

    return field



def get_power_spectra(field1, field2, bins):

    c_l = []

    for zbin1, zbin2 in zip(field1, field2):

        w = nmt.NmtWorkspaceFlat()
        w.compute_coupling_matrix(zbin1, zbin2, bins)

        coupled_c_l = nmt.compute_coupled_cell_flat(zbin1, zbin2)
        c_l.append(w.decouple_cell(coupled_c_l))

    return c_l





def get_power_spectra_all(catfolder):

    nmap, nmap_split1, nmap_split2 = get_nmaps_split(catfolder)

    field0 = [get_nmap(thisnmap) for thisnmap in nmap]
    field1 = [get_nmap(thisnmap) for thisnmap in nmap]
    field2 = [get_nmap(thisnmap) for thisnmap in nmap]

    b = nmt.NmtBinFlat(ell_bins[:-1], ell_bins[1:])

    ells_uncoupled = b.get_effective_ells()

    pair_list = [(field0, field0), (field1, field1), (field1, field2), (field2, field2)]
    
    cl_list = []

    for thisfield1, thisfield2 in pair_list:

        cl_list.append(get_power_spectra(thisfield1, thisfield2, b))

    return ells_uncoupled, cl_list



def plot_power_spectra(catfolder):

    ells_uncoupled, cl_list = get_power_spectra_all(catfolder)

    fig, subplots = plt.subplots(2,2, figsize = (8,8))

    subplots = subplots.flatten()
    names = ['0,0', 's1,s1', 's1,s2', 's2,s2']

    for x, (sp, cl, name) in enumerate(zip(subplots, cl_list, names)):
        sp.plot(ells_uncoupled, cl)
        sp.text(0.98, 0.02, name, transform = sp.transAxes, fontsize = 15)


    












def plot_split_maps(save = False):

    fullmap, splitmap1, splitmap2 = get_nmaps_split()

    fig, ax = plt.subplots(len(newmaps), 2, figsize = (25,30), subplot_kw = {'projection': fsk.wcs})
    for binmap1, binmap2, (axis1, axis2) in zip(newmaps_split1, newmaps_split2, ax):
        fsk.view_map(binmap1, ax = axis1, addColorbar = False)
        fsk.view_map(binmap2, ax = axis2, addColorbar = False)
    
    if save:
        savefig('./nmap_split.pdf', bbox_inches = 'tight')



class Tracer(object) :
    def __init__(self,hdu_list,i_bin,fsk,mask_binary,masked_fraction,contaminants=None) :
        """
        Class used to define the information stored in each of the number density maps generated by CatMapper, which are then transformed into overdensity maps.
        :param hdu_list: list of FITS HDUs containing the number density maps.
        :param i_bin: which redshift bin to consider.
        :param fsk: flatmaps.FlatSkyInfo object defining the geometry of the maps.
        :param mask_binary: binary mask (which pixels to consider and which not to).
        :param masked_fraction: masked fraction map.
        :param contaminants: list of possible contaminant maps to deproject.
        
        This class then stores a number of data objects, the most important one being a pymaster `NmtFieldFlat` ready to use in power spectrum estimation.
        """
        #Read numbers map
        self.fsk,nmap=read_flat_map(None,hdu=hdu_list[2*i_bin])
        compare_infos(fsk,self.fsk)

        #Read N(z)
        self.nz_data=hdu_list[2*i_bin+1].data.copy()

        #Make sure other maps are compatible
        if not self.fsk.is_map_compatible(mask_binary) :
            raise ValueError("Mask size is incompatible")
        if not self.fsk.is_map_compatible(masked_fraction) :
            raise ValueError("Mask size is incompatible")
        if contaminants is not None :
            for ic,c in enumerate(contaminants) :
                if not self.fsk.is_map_compatible(c) :
                    raise ValueError("%d-th contaminant template is incompatible"%ic)
          
        #Translate into delta map
        self.masked_fraction=masked_fraction
        self.weight=masked_fraction*mask_binary
        goodpix=np.where(mask_binary>0.1)[0]
        self.goodpix=goodpix
        self.mask_binary=mask_binary
        self.Ngal = np.sum(nmap*mask_binary)
        ndens=np.sum(nmap*mask_binary)/np.sum(self.weight)
        self.ndens_perad=ndens/(np.radians(self.fsk.dx)*np.radians(self.fsk.dy))
        self.delta=np.zeros_like(self.weight)
        self.delta[goodpix]=nmap[goodpix]/(ndens*masked_fraction[goodpix])-1

        #Reshape contaminants
        conts=None
        if contaminants is not None :
            conts=[[c.reshape([self.fsk.ny,self.fsk.nx])] for c in contaminants]

        #Form NaMaster field
        self.field=nmt.NmtFieldFlat(np.radians(self.fsk.lx),np.radians(self.fsk.ly),
                                    self.weight.reshape([self.fsk.ny,self.fsk.nx]),
                                    [self.delta.reshape([self.fsk.ny,self.fsk.nx])],
templates=conts)