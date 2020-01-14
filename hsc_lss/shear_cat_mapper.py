from ceci import PipelineStage
from .types import FitsFile,ASCIIFile
import numpy as np
from .flatmaps import read_flat_map
from .map_utils import createSpin2Map
from astropy.io import fits

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ShearCatMapper(PipelineStage) :
    name="ShearCatMapper"
    inputs=[('clean_catalog', FitsFile), ('masked_fraction', FitsFile)]
    outputs=[('gamma_maps', FitsFile)]
    config_options={'pz_code':'ephor_ab', 'pz_mark':'best',
                    'pz_bins':[0.15,0.50,0.75,1.00,1.50], 'nz_bin_num':200,
                    'nz_bin_max':3.0, 'shearrot': 'flipuq'}
    
    def get_gamma_maps(self, cat):
        """
        Get gamma1, gamma2 maps and corresponding mask from catalog.
        :param cat:
        :return:
        """

        if not 'ishape_hsm_regauss_e1_calib' in cat.dtype.names:
            raise RuntimeError('get_gamma_maps must be called with calibrated shear catalog. Aborting.')
        maps = []

        for zi, zf in zip(self.zi_arr, self.zf_arr) :
            msk_bin = (cat[self.column_mark]<=zf) & (cat[self.column_mark]>zi)
            subcat = cat[msk_bin]
            gammamaps, gammamasks = createSpin2Map(subcat['ra'], subcat['dec'], subcat['ishape_hsm_regauss_e1_calib'], \
                                     subcat['ishape_hsm_regauss_e2_calib'], self.fsk, \
                                     weights=subcat['ishape_hsm_regauss_derived_shape_weight'], \
                                     shearrot=self.config['shearrot'])
            maps_combined = [gammamaps, gammamasks]
            maps.append(maps_combined)

        return maps

    def parse_input(self):
        """
        Check config parameters for consistency
        """
        #Parse input params
        if self.config['pz_code']=='ephor_ab' :
            self.pz_code='eab'
        elif self.config['pz_code']=='frankenz' :
            self.pz_code='frz'
        elif self.config['pz_code']=='nnpz' :
            self.pz_code='nnz'
        else :
            raise KeyError("Photo-z method "+self.config['pz_code']+
                           " unavailable. Choose ephor_ab, frankenz or nnpz")

        if self.config['pz_mark']  not in ['best','mean','mode','mc'] :
            raise KeyError("Photo-z mark "+self.config['pz_mark']+
                           " unavailable. Choose between best, mean, mode and mc")
        self.column_mark='pz_'+self.config['pz_mark']+'_'+self.pz_code

    def run(self):
        """
        Main routine. This stage:
        - Creates gamma1, gamma2 maps and corresponding masks from the reduced catalog for a set of redshift bins.
        - Stores the above into a single FITS file.
        """

        logger.info("Reading masked fraction from {}.".format(self.get_input("masked_fraction")))
        self.fsk, _ = read_flat_map(self.get_input("masked_fraction"))

        logger.info("Reading calibrated shear catalog from {}.".format(self.get_input('clean_catalog')))
        cat = fits.open(self.get_input('clean_catalog'))[1].data

        logger.info("Creating shear maps and corresponding masks.")
        gammamaps = self.get_gamma_maps(cat)

        print("Writing output to {}.".format(self.get_output('gamma_maps')))
        header = self.fsk.wcs.to_header()
        hdus = []
        for im, m_list in enumerate(gammamaps) :
            # Maps
            if im == 0 :
                head = header.copy()
                head['DESCR'] = ('gamma1, bin %d'%(im+1), 'Description')
                hdu = fits.PrimaryHDU(data=m_list[0][0].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head = header.copy()
                head['DESCR'] = ('gamma2, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[0][1].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head = header.copy()
                head['DESCR'] = ('gamma weight mask, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[1][0].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head['DESCR'] = ('gamma binary mask, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[1][1].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
            else:
                head = header.copy()
                head['DESCR'] = ('gamma1, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[0][0].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head = header.copy()
                head['DESCR'] = ('gamma2, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[0][1].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head = header.copy()
                head['DESCR'] = ('gamma weight mask, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[1][0].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)
                head['DESCR'] = ('gamma binary mask, bin %d'%(im+1), 'Description')
                hdu = fits.ImageHDU(data=m_list[1][1].reshape([self.fsk.ny,self.fsk.nx]), header=head)
                hdus.append(hdu)

        hdulist = fits.HDUList(hdus)
        hdulist.writeto(self.get_output('gamma_maps'), overwrite=True)

if __name__ == '__main__':
    cls = PipelineStage.main()
