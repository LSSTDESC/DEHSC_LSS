import numpy as np
import flatmaps as fm
from shapely.geometry.polygon import Polygon
from shapely.prepared import prep
import sys as sys
import matplotlib.pyplot as plt
from astropy.io import fits

print "Reading map"
fsk,mp=fm.read_flat_map("/global/cscratch1/sd/damonge/HSC/HSC_processed/WIDE_AEGIS/WIDE_AEGIS_MaskedFraction.fits")
print fsk.pix2pos(np.array([0,fsk.nx-1,fsk.nx*(fsk.ny-1),fsk.nx*fsk.ny-1]))

print "Reading frames"
data=fits.open("/global/cscratch1/sd/damonge/HSC/HSC_processed/HSC_WIDE_frames_proc.fits")[1].data

print "Computing frame coords"
ix_ll,iy_ll,in_ll=fsk.pos2pix2d(data['llcra'],data['llcdecl'])
ix_ul,iy_ul,in_ul=fsk.pos2pix2d(data['ulcra'],data['ulcdecl'])
ix_ur,iy_ur,in_ur=fsk.pos2pix2d(data['urcra'],data['urcdecl'])
ix_lr,iy_lr,in_lr=fsk.pos2pix2d(data['lrcra'],data['lrcdecl'])
#Whether they're inside the field at all
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

print "Polygons for map"
polypix=np.array([[Polygon([(ix,iy),(ix,iy+1),(ix+1,iy+1),(ix+1,iy)]) for ix in np.arange(fsk.nx)] for iy in np.arange(fsk.ny)])#.flatten()
indpix=np.arange(fsk.nx*fsk.ny).reshape([fsk.ny,fsk.nx])

print "Getting pixel intersects and areas"
pix_indices=[]
pix_areas=[]
nvisits=np.zeros_like(mp)
for ip,pfr in enumerate(polyframe) :
    print ip
    def get_intersect_area(px) :
        return pfr.intersection(px).area
    c=np.array(pfr.exterior.coords)
    ixmin=max(int(np.amin(c[:,0]))-1,0)
    iymin=max(int(np.amin(c[:,1]))-1,0)
    ixmax=min(int(np.amax(c[:,0]))+1,fsk.nx)
    iymax=min(int(np.amax(c[:,1]))+1,fsk.ny)
    pix_in_range=(polypix[iymin:iymax,:][:,ixmin:ixmax]).flatten()
    ipix_in_range=(indpix[iymin:iymax,:][:,ixmin:ixmax]).flatten()
    pprep=prep(pfr)
    touched=map(pprep.intersects,pix_in_range)
    indices=ipix_in_range[touched]
#    polint=map(pfr.intersection,pix_in_range[touched])
    areas=map(get_intersect_area,pix_in_range[touched])
    print np.amin(areas)
    nvisits[indices]+=1
#    pix_areas.append(areas)
    pix_indices.append(indices)

plt.imshow(nvisits.reshape([fsk.ny,fsk.nx]),extent=(0,fsk.nx,0,fsk.ny),
          origin='lower',interpolation='nearest')
plt.show()
exit(1)
def plot_area(ax,f,col='k') :
    c=np.array(f.exterior.coords)
    ax.fill(c[:,0],c[:,1],col,alpha=0.01)

def plot_path(ax,f,col='k') :
    c=np.array(f.exterior.coords)
    ax.plot(c[:,0],c[:,1],col+'-')

print "Plotting"
plt.figure()
ax=plt.gca()
ax.imshow(mp.reshape([fsk.ny,fsk.nx]),extent=(0,fsk.nx,0,fsk.ny),
          origin='lower',interpolation='nearest')

#plt.figure()
#ax=plt.gca()
for ip,p in enumerate(polyframe[:]) :
    print ip
    plot_area(ax,p,col='b')
#    if np.amin([iy_ll[ip],iy_ul[ip],iy_ur[ip],iy_lr[ip]])<-100 :
#        print np.array(p.exterior.coords),in_ll[ip],in_ul[ip],in_ur[ip],in_lr[ip],is_in[ip]
#        print ' ',data['llcra'][ip],data['llcdecl'][ip]
#        print ' ',data['ulcra'][ip],data['ulcdecl'][ip]
#        print ' ',data['urcra'][ip],data['urcdecl'][ip]
#        print ' ',data['lrcra'][ip],data['lrcdecl'][ip]
#ax.set_xlim([0,fsk.nx])
#ax.set_ylim([0,fsk.ny])

plot_path(ax,polyfield)

plt.show()
