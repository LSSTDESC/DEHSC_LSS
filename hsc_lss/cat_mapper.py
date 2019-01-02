from ceci import PipelineStage
from .types import FitsFile,ASCIIFile
#from astropy.table import Table,vstack
import numpy as np
#from .flatmaps import FlatMapInfo
#from .map_utils import createCountsMap, createMeanStdMaps, createMask, removeDisconnected
#from .estDepth import get_depth
from astropy.io import fits

class CatMapper(PipelineStage) :
    name="CatMapper"
    inputs=[('clean_catalog',FitsFile),('dust_map',FitsFile),('star_map',FitsFile),
            ('masked_fraction',FitsFile),('depth_map',FitsFile),
            ('ccdtemp_maps',FitsFile),('airmass_maps',FitsFile),('exptime_maps',FitsFile),
            ('skylevel_maps',FitsFile),('sigma_sky_maps',FitsFile),('seeing_maps',FitsFile),
            ('ellipt_maps',FitsFile),('nvisit_maps',FitsFile),('cosmos_weights',FitsFile),
            ('pdf_matched',ASCIIFile)]
    outputs=[('ngal_maps',FitsFile)]
    config_options={'mask_type':'sirius','pz_code':'ephor_ab','pz_mark':'best',
                    'pz_bins':[0.15,0.50,0.75,1.00,1.50],'nz_bin_num':200,
                    'nz_bin_max':3.0,}

    def get_nz_cosmos(self) :
        weights_file=fits.open(self.get_input('cosmos_weights'))[1].data

        pzs=[]
        for zi,zf in zip(self.zi_arr,self.zf_arr) :
            msk_cosmos=(weights_file[self.column_mark]<=zf) & (weights_file[self.column_mark]>zi)
            hz,bz=np.histogram(weights_file[msk_cosmos]['PHOTOZ'],
                               bins=self.config['nz_bin_num'],
                               range=[0.,self.config['nz_bin_max']],
                               weights=weights_file[msk_cosmos]['weight'])
            hnz,bnz=np.histogram(weights_file[msk_cosmos]['PHOTOZ'],
                                 bins=self.config['nz_bin_num'],
                                 range=[0.,self.config['nz_bin_max']])
            ehz=np.zeros(len(hnz)); ehz[hnz>0]=(hz[hnz>0]+0.)/np.sqrt(hnz[hnz>0]+0.)
            pzs.append([bz[:-1],bz[1:],(hz+0.)/np.sum(hz+0.),ehz])
        pzs=np.array(pzs)
        return pzs

    def run(self) :
        #Parse input params
        if self.config['pz_code']=='ephor_ab' :
            pz_code='eab'
        elif self.config['pz_code']=='frankenz' :
            pz_code='frz'
        elif self.config['pz_code']=='nnpz' :
            pz_code='nnz'
        else :
            raise KeyError("Photo-z method "+self.config['pz_code']+
                           " unavailable. Choose ephor_ab, frankenz or nnpz")

        if self.config['pz_mark']  not in ['best','mean','mode','mc'] :
            raise KeyError("Photo-z mark "+self.config['pz_mark']+
                           " unavailable. Choose between best, mean, mode and mc")
        self.column_mark='pz_'+self.config['pz_mark']+'_'+pz_code
        print(self.column_mark)

        print("Reading catalog")
        cat=fits.open(self.get_input('clean_catalog'))[1].data
        #Remove masked objects
        if self.config['mask_type']=='arcturus' :
            msk=cat['mask_Arcturus'].astype(bool)
        elif self.config['mask_type']=='sirius' :
            msk=np.logical_not(cat['iflags_pixel_bright_object_center'])
            msk*=np.logical_not(cat['iflags_pixel_bright_object_any'])
        else :
            raise KeyError("Mask type "+self.config['mask_type']+
                           " not supported. Choose arcturus or sirius")
        cat=cat[msk]

        print("Reading pdf filenames")
        data_syst=np.genfromtxt(self.get_input('pdf_matched'),
                                dtype=[('pzname','|U8'),('fname','|U256')])
        self.pdf_files={n:fn for n,fn in zip(data_syst['pzname'],data_syst['fname'])}
        
        print("Parsing photo-z bins")
        self.zi_arr=self.config['pz_bins'][:-1]
        self.zf_arr=self.config['pz_bins'][ 1:]
        self.nbins=len(self.zi_arr)

        print("Getting COSMOS N(z)s")
        pzs_cosmos=self.get_nz_cosmos()
        import matplotlib.pyplot as plt
        for d in pzs_cosmos :
            plt.plot(0.5*(d[0]+d[1]),d[2])
        plt.show()

        print("Getting pdf stacks")
        #for n in pdf_files.keys() :
        #    fn=pdf_files[n]
        #    f=fits.open(fn)
        #    p=f[1].data['pdf'][msk]
        #    z=f[2].data['bins']
        #    print(n,p.shape,z.shape)
        #print(pdf_files)

        exit(1)
        ####

if __name__ == '__main__':
    cls = PipelineStage.main()
