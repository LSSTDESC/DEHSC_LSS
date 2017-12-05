import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

class FlatMapInfo(object) :
    def __init__(self,x_range,y_range,nx=None,ny=None,dx=None,dy=None) :
        """
        Creates a flat map
        x_range : [x_i,x_f] range in the x axis covered by the map
        y_range : [y_i,y_f] range in the y axis covered by the map
        nx,ny : Number of pixels in the x/y axes. If None, dx/dy must be provided
        dx,dy : Resolution in the x/y axes. If None, nx/ny must be provided
        """
        self.x0=x_range[0]
        self.xf=x_range[1]
        self.lx=self.xf-self.x0
        self.y0=y_range[0]
        self.yf=y_range[1]
        self.ly=self.yf-self.y0
        
        if nx is None and dx is None :
            raise ValueError("Must provide either nx or dx")

        if ny is None and dy is None :
            raise ValueError("Must provide either ny or dy")

        if nx is None :
            self.nx=int(self.lx/dx)+1
        else :
            self.nx=nx
        self.dx=self.lx/self.nx

        if ny is None :
            self.ny=int(self.ly/dy)+1
        else :
            self.ny=ny
        self.dy=self.ly/self.ny

        self.npix=self.nx*self.ny
        
    def get_dims(self) :
        """
        Returns map size
        """
        return [self.ny,self.nx]
        
    def get_size(self) :
        """
        Returns map size
        """
        return self.npix

    def pos2pix(self,x,y) :
        """
        Returns pixel indices for arrays of x and y coordinates.
        Will return -1 if (x,y) lies outside the map
        """
        x=np.asarray(x)
        scalar_input=False
        if x.ndim==0 :
            x=x[None]
            scalar_input=True

        y=np.asarray(y)
        if y.ndim==0 :
            y=y[None]

        if len(x)!=len(y) :
            raise ValueError("x and y must have the same size!")

        ix=np.floor((x-self.x0)/self.dx).astype(int)
        ix_out=np.where(np.logical_or(ix<0,ix>=self.nx))[0]

        iy=np.floor((y-self.y0)/self.dy).astype(int)
        iy_out=np.where(np.logical_or(iy<0,iy>=self.ny))[0]

        ipix=ix+self.nx*iy
        ipix[ix_out]=-1
        ipix[iy_out]=-1

        if scalar_input :
            return np.squeeze(ipix)
        return ipix

    def pix2pos(self,ipix) :
        """
        Returns x,y coordinates of pixel centres for a set of pixel indices.
        """
        ipix=np.asarray(ipix)
        scalar_input=False
        if ipix.ndim==0 :
            ipix=ipix[None]
            scalar_input=True

        i_out=np.where(np.logical_or(ipix<0,ipix>=self.npix))[0]
        if len(i_out)>0 :
            print ipix[i_out]
            raise ValueError("Pixels outside of range")

        ix=ipix%self.nx
        ioff=np.array(ipix-ix)
        iy=ioff.astype(int)/(int(self.nx))

        x=self.x0+(ix+0.5)*self.dx
        y=self.y0+(iy+0.5)*self.dy

        if scalar_input :
            return np.squeeze(x),np.squeeze(y)
        return x,y

    def get_empty_map(self) :
        """
        Returns a map full of zeros
        """
        return np.zeros(self.npix,dtype=float)

    def view_map(self,map_in,ax=None,xlabel='x',ylabel='y',title=None) :
        """
        Plots a 2D map (passed as a flattened array)
        """
        if len(map_in)!=self.npix :
            raise ValueError("Input map doesn't have the correct size")

	# set up the colorbar
	cmap = cm.magma
	cmap.set_under("w")
	# min, max
	median= np.median(map_in)
        stddev= np.std(map_in)
        colorMin= median-1.5*stddev
        colorMax= median+1.5*stddev

        if ax is None :
            plt.figure()
            ax=plt.gca()
        if title is not None :
            ax.set_title(title,fontsize=15)
        image= ax.imshow(map_in.reshape([self.ny,self.nx]),
			origin='lower', interpolation='nearest',
			aspect='auto', extent=[self.x0,self.xf,self.y0,self.yf],
			vmin= colorMin, vmax= colorMax, cmap= cmap)
	plt.colorbar(image)
        ax.set_xlabel(xlabel,fontsize=15)
        ax.set_ylabel(ylabel,fontsize=15)

    def write_flat_map(self,filename,maps) :
        """
        Saves a set of maps in npz format.
        We'll try to implement other more standard formats with proper WCS coordinates etc. ASAP.
        """
        if maps.ndim<1 :
            raise ValueError("Must supply at least one map")
        if maps.ndim==1 :
            maps=np.array([maps])
        if len(maps[0])!=self.npix :
            raise ValueError("Map doesn't conform to this pixelization")
            
        np.savez(filename,x_range=[self.x0,self.xf],y_range=[self.y0,self.yf],nx=self.nx,ny=self.ny,
                 maps=maps)



def read_flat_map(filename,i_map=0) :
    """
    Reads a flat-sky map and the details of its pixelization scheme.
    The latter are returned as a FlatMapInfo object.
    i_map : map to read. If -1, all maps will be read.
    """
    data=np.load(filename)

    fmi=FlatMapInfo(data['x_range'],data['y_range'],nx=data['nx'],ny=data['ny'])
    if i_map==-1 :
        i_map=np.arange(len(data['maps']))
    return fmi,data['maps'][i_map]
