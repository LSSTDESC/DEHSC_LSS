from ceci import PipelineStage
from .types import FitsFile,ASCIIFile
import numpy as np
from .flatmaps import read_flat_map
from .map_utils import createCountsMap
from astropy.io import fits

class CatMapper(PipelineStage) :
    name="CatMapper"
    inputs=[('clean_catalog',FitsFile),('masked_fraction',FitsFile),
            ('cosmos_weights',FitsFile),('pdf_matched',ASCIIFile)]
    outputs=[('ngal_maps',FitsFile)]
    config_options={'mask_type':'sirius','pz_code':'ephor_ab','pz_mark':'best',
                    'pz_bins':[0.15,0.50,0.75,1.00,1.50],'nz_bin_num':200,
                    'nz_bin_max':3.0,}
    
    def get_nmaps(self,cat) :
        maps=[]

        for zi,zf in zip(self.zi_arr,self.zf_arr) :
            msk_bin=(cat[self.column_mark]<=zf) & (cat[self.column_mark]>zi)
            subcat=cat[msk_bin]
            nmap=createCountsMap(subcat['ra'],subcat['dec'],self.fsk)
            maps.append(nmap)
        return np.array(maps)

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
        return np.array(pzs)

    def get_nz_stack(self,cat,codename) :
        from scipy.interpolate import interp1d

        f=fits.open(self.pdf_files[codename])
        p=f[1].data['pdf'][self.msk]
        z=f[2].data['bins']

        z_all=np.linspace(0.,self.config['nz_bin_max'],self.config['nz_bin_num']+1)
        z0=z_all[:-1]; z1=z_all[1:]; zm=0.5*(z0+z1)
        pzs=[]
        for zi,zf in zip(self.zi_arr,self.zf_arr) :
            msk_bin=(cat[self.column_mark]<=zf) & (cat[self.column_mark]>zi)
            hz_orig=np.sum(p[msk_bin],axis=0)
            hz_orig/=np.sum(hz_orig)
            hzf=interp1d(z,hz_orig,bounds_error=False,fill_value=0.)
            hzm=hzf(zm);
            
            pzs.append([z0,z1,hzm/np.sum(hzm)])
        return np.array(pzs)
            
    def parse_input(self) :
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

    def run(self) :
        self.parse_input()
        
        print("Reading masked fraction")
        self.fsk,_=read_flat_map(self.get_input("masked_fraction"))

        print("Reading catalog")
        cat=fits.open(self.get_input('clean_catalog'))[1].data
        #Remove masked objects
        if self.config['mask_type']=='arcturus' :
            self.msk=cat['mask_Arcturus'].astype(bool)
        elif self.config['mask_type']=='sirius' :
            self.msk=np.logical_not(cat['iflags_pixel_bright_object_center'])
            self.msk*=np.logical_not(cat['iflags_pixel_bright_object_any'])
        else :
            raise KeyError("Mask type "+self.config['mask_type']+
                           " not supported. Choose arcturus or sirius")
        cat=cat[self.msk]

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

        print("Getting pdf stacks")
        pzs_stack={}
        for n in self.pdf_files.keys() :
            pzs_stack[n]=self.get_nz_stack(cat,n)

        print("Getting number count maps")
        n_maps=self.get_nmaps(cat)

        print("Writing output")
        header=self.fsk.wcs.to_header()
        hdus=[]
        for im,m in enumerate(n_maps) :
            #Map
            head=header.copy()
            head['DESCR']=('Ngal, bin %d'%(im+1),'Description')
            if im==0 :
                hdu=fits.PrimaryHDU(data=m.reshape([self.fsk.ny,self.fsk.nx]),header=head)
            else :
                hdu=fits.ImageHDU(data=m.reshape([self.fsk.ny,self.fsk.nx]),header=head)
            hdus.append(hdu)
            
            #Nz
            cols=[fits.Column(name='z_i',array=pzs_cosmos[im,0,:],format='E'),
                  fits.Column(name='z_f',array=pzs_cosmos[im,1,:],format='E'),
                  fits.Column(name='nz_cosmos',array=pzs_cosmos[im,2,:],format='E'),
                  fits.Column(name='enz_cosmos',array=pzs_cosmos[im,3,:],format='E')]
            for n in self.pdf_files.keys() :
                cols.append(fits.Column(name='nz_'+n,array=pzs_stack[n][im,2,:],format='E'))
            hdus.append(fits.BinTableHDU.from_columns(cols))
        hdulist=fits.HDUList(hdus)
        hdulist.writeto(self.get_output('ngal_maps'),overwrite=True)

if __name__ == '__main__':
    cls = PipelineStage.main()
