import numpy as np

import copy
from astropy.table import Table, vstack
from astropy.io import fits
from ceci import PipelineStage
from .types import FitsFile, ASCIIFile

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ShearCalib(PipelineStage):

    name = "ReduceCat"
    inputs = [('raw_data', None)]
    outputs = [('calib_catalog', FitsFile),('R', ASCIIFile),('mhat',ASCIIFile)]
    config_options = {'photoz_method': 'ephor_ab_photoz_best', 'photoz_min':0.3, 'photoz_max': 1.5}

    def _responsivity(self, cat):
        """
        Compute shear responsivity.
        For HSC (see Mandelbaum et al., 2018, arXiv:1705.06745):
        R = 1 - < e_rms^2 >w (Eq. (A1) in Mandelbaum et al., 2018)
        :param cat:
        :return:
        """

        R = 1. - np.average(cat['ishape_hsm_regauss_derived_rms_e']**2, weights=cat['ishape_hsm_regauss_derived_shape_weight'])

        return R

    def _mhat(self, cat):
        """
        Compute multiplicative bias.
        For HSC (see Mandelbaum et al., 2018, arXiv:1705.06745):
        mhat = < m >w (Eq. (A2) in Mandelbaum et al., 2018)
        :param cat:
        :return:
        """

        mhat = np.average(cat['ishape_hsm_regauss_derived_shear_bias_m'], weights=cat['ishape_hsm_regauss_derived_shape_weight'])

        return mhat

    def calibrated_catalog(self, cat, R=None, mhat=None):
        """
        Calibrate shear catalog and add calibrated shear columns to existing catalog.
        For HSC (see Mandelbaum et al., 2018, arXiv:1705.06745):
        gi = 1/(1 + mhat)[ei/(2R) - ci] (Eq. (A6) in Mandelbaum et al., 2018)
        R = 1 - < e_rms^2 >w (Eq. (A1) in Mandelbaum et al., 2018)
        mhat = < m >w (Eq. (A2) in Mandelbaum et al., 2018)
        :param cat:
        :param R:
        :param mhat:
        :return:
        """

        logger.info('Computing calibrated shear catalog.')

        cat_calib = copy.deepcopy(cat)

        if R is None and mhat is None:
            logger.info('Computing R and mhat.')
            R = self._responsivity(cat_calib)
            mhat = self._mhat(cat_calib)

        else:
            logger.info('R and mhat provided.')

        logger.info('R = {}, mhat = {}.'.format(R, mhat))

        e1_corr = 1./(1. + mhat)*(cat_calib['ishape_hsm_regauss_e1']/(2.*R) - cat_calib['ishape_hsm_regauss_derived_shear_bias_c1'])
        e2_corr = 1./(1. + mhat)*(cat_calib['ishape_hsm_regauss_e2']/(2.*R) - cat_calib['ishape_hsm_regauss_derived_shear_bias_c2'])

        # Add these two columns to catalog
        cat_calib = Table(cat_calib)
        cat_calib['ishape_hsm_regauss_e1_calib'] = e1_corr
        cat_calib['ishape_hsm_regauss_e2_calib'] = e2_corr
        cat_calib = np.array(cat_calib)

        logger.info('Columns ishape_hsm_regauss_e1_calib, ishape_hsm_regauss_e2_calib added to shear catalog.')

        return cat_calib, e1_corr, e2_corr

    def pz_cut(self, cat):
        """
        Apply pz cut to catalog.
        :param cat:
        :return:
        """

        logger.info('Applying pz cut to catalog. Using {} with zmin = {}, zmax = {}.'.\
                    format(self.config['photoz_method'], self.config['photoz_min'], self.config['photoz_max']))

        photozmask = (cat[self.config['photoz_method']]>=self.config['photoz_min'])\
                     &(cat[self.config['photoz_method']]<self.config['photoz_max'])

        cat = copy.deepcopy(cat)
        cat = cat[photozmask]

        return cat

    def run(self) :
        """
        Main function.
        This stage:
        - Reduces shear catalog.
        - Computes responsivity and multiplicative bias.
        - Computes calibrated shear components g1, g2 and writes calibrated catalog.
        """

        # Read list of files
        f = open(self.get_input('raw_data'))
        files = [s.strip() for s in f.readlines()]
        f.close()

        #Read catalog
        cat = Table.read(files[0])
        if len(cat) > 1 :
            for fname in files[1:]:
                c = Table.read(fname)
                cat = vstack([cat, c], join_type='exact')

        logger.info('Initial catalog size: {}.'.format(len(cat)))
        cat = self.pz_cut(cat)
        logger.info('Catalog size after pz cut: {}.'.format(len(cat)))

        cat_calib, R, mhat = self.calibrated_catalog(cat)

        ####
        # Write calibrated catalog
        # 1- header
        logger.info('Writing calibrated shear catalog.')
        hdr = fits.Header()
        hdr['R'] = R
        hdr['mhat'] = mhat
        prm_hdu = fits.PrimaryHDU(header=hdr)
        # 2- Catalog
        cat_hdu = fits.table_to_hdu(cat_calib)
        # 3- Actual writing
        hdul = fits.HDUList([prm_hdu,cat_hdu])
        hdul.writeto(self.get_output('clean_catalog'), overwrite=True)

        ####

if __name__ == '__main__':
    cls = PipelineStage.main()