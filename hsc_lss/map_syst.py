from ceci import PipelineStage
from .types import FitsFile
import numpy as np
from .flatmaps import FlatMapInfo, read_flat_map
from .obscond import ObsCond
#from .map_utils import createCountsMap, createMeanStdMaps, createMask, removeDisconnected
#from .estDepth import get_depth
from astropy.io import fits
from shapely.geometry.polygon import Polygon
from shapely.prepared import prep

class MapSyst(PipelineStage) :
    name="MapSyst"
    inputs=[('frames_data',FitsFile),('masked_fraction',FitsFile)]
    outputs=[('ccdtemp_maps',FitsFile),('airmass_maps',FitsFile),('exptime_maps',FitsFile),
             ('skylevel_maps',FitsFile),('sigma_sky_maps',FitsFile),('seeing_maps',FitsFile),
             ('ellipt_maps',FitsFile),('nvisit_maps',FitsFile)]
    config_options={'hist_nbins':40}

    def run(self) :
        bands=['g','r','i','z','y']
        quants=['ccdtemp','airmass','exptime','skylevel','sigma_sky','seeing','ellipt']
        quants_islog=[False,False,False,True,True,False,False]

        print("Reading sample map")
        fsk,mp=read_flat_map(self.get_input('masked_fraction'))

        print("Reading metadata")
        data=fits.open(self.get_input('frames_data'))[1].data

        print("Computing frame coords")
        ix_ll,iy_ll,in_ll=fsk.pos2pix2d(data['llcra'],data['llcdecl'])
        ix_ul,iy_ul,in_ul=fsk.pos2pix2d(data['ulcra'],data['ulcdecl'])
        ix_ur,iy_ur,in_ur=fsk.pos2pix2d(data['urcra'],data['urcdecl'])
        ix_lr,iy_lr,in_lr=fsk.pos2pix2d(data['lrcra'],data['lrcdecl'])
        #Keep only frames that fit inside the field
        is_in=np.logical_or(in_ll,np.logical_or(in_ul,np.logical_or(in_ur,in_lr)))
        data=data[is_in]
        nframes=len(data)
        ix_ll=ix_ll[is_in]; iy_ll=iy_ll[is_in]; 
        ix_ul=ix_ul[is_in]; iy_ul=iy_ul[is_in]; 
        ix_ur=ix_ur[is_in]; iy_ur=iy_ur[is_in]; 
        ix_lr=ix_lr[is_in]; iy_lr=iy_lr[is_in];
        
        print("Building poliygons")
        polyfield=Polygon([(0,0),(0,fsk.ny),(fsk.nx,fsk.ny),(fsk.nx,0)])
        polyframe=np.array([Polygon([(ix_ll[i],iy_ll[i]),
                                     (ix_ul[i],iy_ul[i]),
                                     (ix_ur[i],iy_ur[i]),
                                     (ix_lr[i],iy_lr[i])]) for i in np.arange(len(data))])
        polypix=np.array([[Polygon([(ix,iy),(ix,iy+1),(ix+1,iy+1),(ix+1,iy)]) 
                           for ix in np.arange(fsk.nx)] 
                          for iy in np.arange(fsk.ny)])
        indpix=np.arange(fsk.nx*fsk.ny).reshape([fsk.ny,fsk.nx])

        print("Getting pixel intersects and areas")
        pix_indices=[]
        pix_areas=[]
        percent_next=0
        for ip,pfr in enumerate(polyframe) :
            percent_done=int(100*(ip+0.)/nframes)
            if percent_done==percent_next :
                print("%d%% done"%percent_done)
                percent_next+=10
            c=np.array(pfr.exterior.coords)
            
            #Pick only pixels in range
            ixmin=max(int(np.amin(c[:,0]))-1,0)
            iymin=max(int(np.amin(c[:,1]))-1,0)
            ixmax=min(int(np.amax(c[:,0]))+1,fsk.nx)
            iymax=min(int(np.amax(c[:,1]))+1,fsk.ny)
            pix_in_range=(polypix[iymin:iymax,:][:,ixmin:ixmax]).flatten()
            ipix_in_range=(indpix[iymin:iymax,:][:,ixmin:ixmax]).flatten()

            #Estimate which ones actually have been touched
            pprep=prep(pfr)
            touched=list(map(pprep.intersects,pix_in_range))
            indices=ipix_in_range[touched]
            pix_indices.append(indices)

            #Estimate intersection area
            def get_intersect_area(px) :
                return pfr.intersection(px).area
            areas=list(map(get_intersect_area,pix_in_range[touched]))
            pix_areas.append(areas)

        print("Computing systematics maps")
        #Initialize maps
        nvisits={b:np.zeros_like(mp) for b in bands}
        oc_maps={}
        for l,q in zip(quants_islog,quants) :
            oc_maps[q]={b:ObsCond(q,data[q],fsk.nx,fsk.ny,nbins=self.config['hist_nbins'],is_log=l) 
                        for b in bands}
        #Fill maps    
        for ip,pfr in enumerate(polyframe) :
            band=data['filter'][ip]
            indices=pix_indices[ip]
            areas=pix_areas[ip]
            nvisits[band][pix_indices[ip]]+=pix_areas[ip]
            for quant in quants :
                d=data[quant][ip]
                if d<=-9999. :
                    continue
                ib=oc_maps[quant][band].get_bin_number(d)
                oc_maps[quant][band].map[indices,ib]+=areas

        print("Saving maps")
        #Nvisits
        maps_save=np.array([nvisits[b] for b in bands])
        descripts=np.array(['Nvisits-'+b for b in bands])
        fsk.write_flat_map(self.get_output('nvisit_maps'),maps_save,descripts)
        #Observing conditions
        for q in quants :
            maps_save=np.array([oc_maps[q][b].collapse_map_mean() for b in bands])
            descripts=np.array(['mean '+q+'-'+b for b in bands])
            fsk.write_flat_map(self.get_output(q+'_maps'),maps_save,descripts)


if __name__ == '__main__':
    cls = PipelineStage.main()
