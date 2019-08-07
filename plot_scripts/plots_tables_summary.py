import numpy as np
from astropy.io import fits

# Field summary
def get_field_stats(fieldname):
    cat=fits.open("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+fieldname+"_sirius_i24p5_out/clean_catalog.fits")[1].data
    mskd=fits.open("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+fieldname+"_sirius_i24p5_out/masked_fraction.fits")[0]

    area_pix_deg2=np.fabs(mskd.header['CDELT1']*mskd.header['CDELT2'])
    area_pix_sr=area_pix_deg2*(np.pi/180)**2
    area_deg2=np.sum(mskd.data)*area_pix_deg2

    ng=len(cat)
    print(fieldname+" %d gals, %.1lf deg^2, fsky=%.1lE"%(ng,area_deg2,area_deg2/(4*np.pi*(180/np.pi)**2)))
    return ng,area_deg2
ngt=0
area_deg2t=0
for field in ['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']:
    n,a=get_field_stats(field)
    ngt+=n
    area_deg2t+=a
print("Total: %d gals, %.1lf deg^2, fsky=%.1lE"%(ngt,area_deg2t,area_deg2t/(4*np.pi*(180/np.pi)**2)))

# Bin summary
pz_bins=[0.15,0.50,0.75,1.00,1.50]
ng_arr=[]
nz_arr=[]
for fieldname in ['GAMA09H','GAMA15H','HECTOMAP','VVDS','WIDE12H','XMMLSS']:
    mskd=fits.open("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+fieldname+"_sirius_i24p5_out/masked_fraction.fits")[0].data
    mskbin=np.zeros_like(mskd); mskbin[mskd>0.5]=1
    fng=fits.open("/global/cscratch1/sd/damonge/HSC_ceci/WIDE_"+fieldname+"_sirius_i24p5_out/ngal_maps.fits")
    zs=0.5*(fng[1].data['z_f']+fng[1].data['z_i'])
    ngs=[]
    nzs=[]
    for i in range(4):
        #print(fng[0].data.shape)
        #print(fng[1].data.names)
        ngs.append(np.sum(fng[2*i].data*mskbin))
        nzs.append(fng[2*i+1].data['nz_cosmos'])
    ng_arr.append(ngs)
    nz_arr.append(nzs)
ng_arr=np.array(ng_arr)
nz_arr=np.array(nz_arr)
ng=np.sum(ng_arr,axis=0)
nz=np.sum(ng_arr[:,:,None]*nz_arr,axis=0)
print("z_mean : ")
print(np.sum(nz[:,:]*zs[None,:],axis=1)/np.sum(nz,axis=1),np.sum(ng[:,None]*nz*zs[None,:])/np.sum(ng[:,None]*nz))
print("n_gal : ")
print(ng,np.sum(ng))
