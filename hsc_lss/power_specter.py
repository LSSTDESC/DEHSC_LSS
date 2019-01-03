from ceci import PipelineStage
from .types import FitsFile,ASCIIFile
import numpy as np
from .flatmaps import read_flat_map,compare_infos
from astropy.io import fits

class PowerSpecter(PipelineStage) :
    name="PowerSpecter"
    inputs=[('masked_fraction',FitsFile),('ngal_maps',FitsFile),
            ('dust_map',FitsFile),('star_map',FitsFile),('depth_map',FitsFile),
            ('ccdtemp_maps',FitsFile),('airmass_maps',FitsFile),('exptime_maps',FitsFile),
            ('skylevel_maps',FitsFile),('sigma_sky_maps',FitsFile),('seeing_maps',FitsFile),
            ('ellipt_maps',FitsFile),('nvisit_maps',FitsFile),('cosmos_weights',FitsFile)]
    outputs=[('power_spectra_wdpj',FitsFile),('power_spectra_wodpj',FitsFile)]
    config_options={'ell_bpws':[100.0,200.0,300.0,
                                400.0,600.0,800.0,
                                1000.0,1400.0,1800.0,
                                2200.0,3000.0,3800.0,
                                4600.0,6200.0,7800.0,
                                9400.0,12600.0,12600.0,15800.0],
                    'mcm_file':'NONE','mcm_covar_file':'NONE','depth_cut':24.5,
                    'mask_thr':0.5,'guess_spectrum':'NONE','gaus_covar_type':'analytic',
                    'add_ssc':False,'mask_systematics':False,'noise_bias_type':'analytic'}
            
    def parse_input(self) :
        #Parse input params

        return

    def run(self) :
        self.parse_input()

        import matplotlib.pyplot as plt

        print("Reading mask")
        #Depth-based mask
        self.fsk,mp_depth=read_flat_map(self.get_input("depth_map"),i_map=0)
        mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
        msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=self.config['depth_cut']]=1
        self.fsk.view_map(msk_depth);

        fskb,mskfrac=read_flat_map(self.get_input("masked_fraction"),i_map=0)
        compare_infos(self.fsk,fskb)
        self.fsk.view_map(mskfrac)
        
        #Create binary mask (fraction>threshold and depth req.)
        msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>self.config['mask_thr']]=1
        msk_bi=msk_bo*msk_depth

        self.fsk.view_map(msk_bi*mskfrac)

        if self.config['mask_systematics'] :
            raise NotImplementedError("TODO: implement systematics masking")

        plt.show(); exit(1)

if __name__ == '__main__':
    cls = PipelineStage.main()
