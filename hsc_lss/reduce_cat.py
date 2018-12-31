from ceci import PipelineStage
from .types import FitsFile
from astropy.table import Table,vstack
import numpy as np
from .flatmaps import FlatMapInfo
from .map_utils import createCountsMap, createMeanStdMaps, createMask, removeDisconnected

class ReduceCat(PipelineStage) :
    name="ReduceCat"
    inputs=[('raw_data',None)]
    outputs=[('clean_catalog',FitsFile),('dust_map',FitsFile),('star_map',FitsFile),
             ('bo_mask',FitsFile),('masked_fraction',FitsFile),('depth_map',FitsFile)]
    config_options={'min_snr':10.,'depth_cut':24.5,'res':0.0285,
                    'res_bo':0.003,'pad':0.1,'band':'i','depth_method':2,
                    'flat_project':'CAR','mask_type':'sirius'}
    bands=['g','r','i','z','y']

    def make_dust_map(self,cat,fsk) :
        print("Creating dust map")
        dustmaps=[]
        dustdesc=[]
        for b in self.bands :
            m,s=createMeanStdMaps(cat['ra'],cat['dec'],cat['a_'+b],fsk)
            dustmaps.append(m)
            dustdesc.append('Dust, '+b+'-band')
        return dustmaps,dustdesc

    def make_star_map(self,cat,fsk,sel) :
        print("Creating star map")
        mstar=createCountsMap(cat['ra'][sel],cat['dec'][sel],fsk)+0.
        descstar='Stars, '+self.config['band']+'<%.2lf'%(self.config['depth_cut'])
        return mstar,descstar

    def make_bo_mask(self,cat,fsk) :
        print("Generating bright-object mask")
        if self.config['mask_type']=='arcturus' :
            flags_mask=[~cat['mask_Arcturus'].astype(bool)]
        elif self.config['mask_type']=='sirius' :
            flags_mask=[cat['iflags_pixel_bright_object_center'],
                        cat['iflags_pixel_bright_object_any']]
        else :
            raise ValueError('Mask type '+self.config['mask_type']+' not supported')
        mask_bo,fsg=createMask(cat['ra'],cat['dec'],flags_mask,fsk,self.config['res_bo'])
        return mask_bo,fsg

    def make_masked_fraction(self,cat,fsk) :
        print("Generating masked fraction map")
        masked=np.ones(len(cat))
        if self.config['mask_type']=='arcturus' :
            masked*=cat['mask_Arcturus']
        elif self.config['mask_type']=='sirius' :
            masked*=np.logical_not(cat['iflags_pixel_bright_object_center'])
            masked*=np.logical_not(cat['iflags_pixel_bright_object_any'])
        else :
            raise ValueError('Mask type '+self.config['mask_type']+' not supported')
        masked_fraction,_=createMeanStdMaps(cat['ra'],cat['dec'],masked,fsk)
        masked_fraction_cont=removeDisconnected(masked_fraction,fsk)
        return masked_fraction_cont

    def run(self) :
        #Read list of files
        f=open(self.get_input('raw_data'))
        files=[s.strip() for s in f.readlines()]
        f.close()

        #Read catalog
        cat=Table.read(files[0])
        if len(cat)>1 :
            for fname in files[1:] :
                c=Table.read(fname)
                cat=vstack([cat,c],join_type='exact')

        if self.config['band'] not in self.bands :
            raise ValueError("Band "+self.config['band']+" not available")
            
        # Clean nulls and nans
        print("Basic cleanup")
        sel=np.ones(len(cat),dtype=bool)
        names=[n for n in cat.keys()]
        isnull_names=[]
        for key in cat.keys() :
            if key.__contains__('isnull') :
                sel[cat[key]]=0
                isnull_names.append(key)
            else :
                if not key.startswith("pz_") : #Keep photo-z's even if they're NaNs
                    sel[np.isnan(cat[key])]=0
        print("Will drop %d rows"%(len(sel)-np.sum(sel)))
        cat.remove_columns(isnull_names)
        cat.remove_rows(~sel)

        fsk=FlatMapInfo.from_coords(cat['ra'],cat['dec'],self.config['res'],
                                    pad=self.config['pad']/self.config['res'],
                                    projection=self.config['flat_project'])

        #Collect sample cuts
        sel_maglim=np.ones(len(cat),dtype=bool);
        sel_maglim[cat['%scmodel_mag'%self.config['band']]-
                   cat['a_%s'%self.config['band']]>self.config['depth_cut']]=0
        # Blending
        sel_blended=np.ones(len(cat),dtype=bool);
        sel_blended[cat['iblendedness_abs_flux']>=0.42169650342]=0 #abs_flux<10^-0.375
        # S/N in i
        sel_fluxcut_i=np.ones(len(cat),dtype=bool);
        sel_fluxcut_i[cat['icmodel_flux']<10*cat['icmodel_flux_err']]=0
        # S/N in g
        sel_fluxcut_g=np.ones(len(cat),dtype=int);
        #TODO: this is a bug leftover from the first pipeline
        sel_fluxcut_i[cat['gcmodel_flux']<5*cat['gcmodel_flux_err']]=0
        # S/N in r
        sel_fluxcut_r=np.ones(len(cat),dtype=int);
        sel_fluxcut_i[cat['rcmodel_flux']<5*cat['rcmodel_flux_err']]=0
        # S/N in z
        sel_fluxcut_z=np.ones(len(cat),dtype=int);
        sel_fluxcut_i[cat['zcmodel_flux']<5*cat['zcmodel_flux_err']]=0
        # S/N in y
        sel_fluxcut_y=np.ones(len(cat),dtype=int);
        sel_fluxcut_i[cat['ycmodel_flux']<5*cat['ycmodel_flux_err']]=0
        # S/N in grzy (at least 2 pass)
        sel_fluxcut_grzy=(sel_fluxcut_g+sel_fluxcut_r+sel_fluxcut_z+sel_fluxcut_y>=2)
        # Overall S/N
        sel_fluxcut=sel_fluxcut_i*sel_fluxcut_grzy
        # Stars
        sel_stars=np.ones(len(cat),dtype=bool);
        sel_stars[cat['iclassification_extendedness']>0.99]=0
        # Galaxies
        sel_gals =np.ones(len(cat),dtype=bool);
        sel_gals[cat['iclassification_extendedness']<0.99]=0

        ####
        # Generate systematics maps
        # 1- Dust
        dustmaps,dustdesc=self.make_dust_map(cat,fsk)
        fsk.write_flat_map(self.get_output('dust_map'),np.array(dustmaps),descript=dustdesc)

        # 2- Nstar
        #    This needs to be done for stars passing the same cuts as the sample 
        #    (except for the s/g separator)
        # Above magnitude limit
        mstar,descstar=self.make_star_map(cat,fsk,sel_maglim*sel_stars*sel_fluxcut*sel_blended)
        fsk.write_flat_map(self.get_output('star_map'),mstar,descript=descstar)

        
        #Binary BO mask
        mask_bo,fsg=self.make_bo_mask(cat,fsk)
        fsg.write_flat_map(self.get_output('bo_mask'),mask_bo,descript='Bright-object mask')

        #Masked fraction
        masked_fraction_cont=self.make_masked_fraction(cat,fsk)
        fsk.write_flat_map(self.get_output('masked_fraction'),masked_fraction_cont,
                           descript='Masked fraction')


if __name__ == '__main__':
    cls = PipelineStage.main()
