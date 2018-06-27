import numpy as np
import flatmaps as fm
from shapely.geometry.polygon import Polygon
from shapely.prepared import prep
import sys as sys
import matplotlib.pyplot as plt
from astropy.io import fits
from obscond import ObsCond

fname_map_sample="/global/cscratch1/sd/damonge/HSC/HSC_processed/WIDE_AEGIS/WIDE_AEGIS_MaskedFraction.fits"
fname_frames="/global/cscratch1/sd/damonge/HSC/HSC_processed/HSC_WIDE_frames_proc.fits"
output_prefix="/global/cscratch1/sd/damonge/HSC/HSC_processed/WIDE_AEGIS/WIDE_AEGIS"
quants=['ccdtemp','airmass','exptime','skylevel','sigma_sky','seeing','ellipt']
quants_islog=[False,False,False,True,True,False,False]
bands=['g','r','i','z','y']

print "Reading map"
fsk,mp=fm.read_flat_map(fname_map_sample)

print "Reading frames"
data=fits.open(fname_frames)[1].data
def plotstats(q,islog=False) :
    print q,np.amin(data[q][data[q]>-999.]),np.amax(data[q]),np.mean(data[q])
    dplot=data[q][data[q]>-9999.]
    if islog :
        dplot=np.log10(dplot)
    plt.figure(); plt.title(q); plt.hist(dplot,bins=100)

print "Computing frame coords"
ix_ll,iy_ll,in_ll=fsk.pos2pix2d(data['llcra'],data['llcdecl'])
ix_ul,iy_ul,in_ul=fsk.pos2pix2d(data['ulcra'],data['ulcdecl'])
ix_ur,iy_ur,in_ur=fsk.pos2pix2d(data['urcra'],data['urcdecl'])
ix_lr,iy_lr,in_lr=fsk.pos2pix2d(data['lrcra'],data['lrcdecl'])
#Keep only frames that fit inside the field
is_in=np.logical_or(in_ll,np.logical_or(in_ul,np.logical_or(in_ur,in_lr)))
data=data[is_in]
ix_ll=ix_ll[is_in]; iy_ll=iy_ll[is_in]; 
ix_ul=ix_ul[is_in]; iy_ul=iy_ul[is_in]; 
ix_ur=ix_ur[is_in]; iy_ur=iy_ur[is_in]; 
ix_lr=ix_lr[is_in]; iy_lr=iy_lr[is_in];

print "Polygon for field"
polyfield=Polygon([(0,0),(0,fsk.ny),(fsk.nx,fsk.ny),(fsk.nx,0)])

print "Polygons for frames"
polyframe=np.array([Polygon([(ix_ll[i],iy_ll[i]),
                             (ix_ul[i],iy_ul[i]),
                             (ix_ur[i],iy_ur[i]),
                             (ix_lr[i],iy_lr[i])]) for i in np.arange(len(data))])

print "Polygons for pixels"
polypix=np.array([[Polygon([(ix,iy),(ix,iy+1),(ix+1,iy+1),(ix+1,iy)]) 
                   for ix in np.arange(fsk.nx)] 
                  for iy in np.arange(fsk.ny)])
indpix=np.arange(fsk.nx*fsk.ny).reshape([fsk.ny,fsk.nx])

print "Getting pixel intersects and areas"
pix_indices=[]
pix_areas=[]
for ip,pfr in enumerate(polyframe) :
    if ip%100==0 :
        print ip
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
    touched=map(pprep.intersects,pix_in_range)
    indices=ipix_in_range[touched]
    pix_indices.append(indices)

    #Estimate intersection area
    def get_intersect_area(px) :
        return pfr.intersection(px).area
    areas=map(get_intersect_area,pix_in_range[touched])
    pix_areas.append(areas)

print "Initializing systematics maps"
nvisits={b:np.zeros_like(mp) for b in bands}
oc_maps={}
for l,q in zip(quants_islog,quants) :
    oc_maps[q]={b:ObsCond(q,data[q],fsk.nx,fsk.ny,is_log=l) for b in bands}

print "Computing systematics maps"
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

print "Saving maps"
#Nvisits
maps_save=np.array([nvisits[b] for b in bands])
descripts=np.array(['Nvisits-'+b for b in bands])
fsk.write_flat_map(output_prefix+'_oc_nvisit.fits',maps_save,descripts)
#Observing conditions
for q in quants :
    maps_save=np.array([oc_maps[q][b].collapse_map_mean() for b in bands])
    descripts=np.array(['mean '+q+'-'+b for b in bands])
    fsk.write_flat_map(output_prefix+'_oc_'+q+'.fits',maps_save,descripts)
