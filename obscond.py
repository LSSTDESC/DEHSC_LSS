from __future__ import print_funciton
import numpy as np

class ObsCond(object):
    def __init__(self,name,data,nx,ny,nbins=40,is_log=False,cutoff=-9999.) :
        self.name=name
        self.is_log=is_log

        d_use=data[data>cutoff]
        if self.is_log :
            d_use=np.log10(d_use)
        
        self.nbins=nbins
        self.dmin=np.amin(d_use)
        self.dmax=np.amax(d_use)
        #Extra space at the end
        self.dmax+=0.001*(self.dmax-self.dmin)
        self.dd=(self.dmax-self.dmin)/self.nbins
        self.idd=1./self.dd
        self.xarr=self.dmin+(np.arange(self.nbins)+0.5)*self.dd
 
        self.map=np.zeros([ny*nx,nbins])

    def get_bin_number(self,d) :
        if self.is_log :
            x=np.log10(d)
        else :
            x=d

        return int(self.idd*(x-self.dmin))

    def collapse_map_mean(self) :
        norm=np.sum(self.map,axis=1)

        map_out=np.zeros_like(norm)-9999.
        
        #Check for empty pixels
        indgood=np.where(norm>0.)[0]
        map_out[indgood]=np.sum(self.map[indgood]*self.xarr[None,:],axis=1)/norm[indgood]

        if self.is_log :
            map_out[indgood]=10.**map_out[indgood]

        return map_out
