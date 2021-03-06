from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from createMaps import createCountsMap
from optparse import OptionParser
import flatmaps as fm
import pymaster as nmt
import sacc
import sys
import time
import os
from scipy.interpolate import interp1d

def opt_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
#Options
parser.add_option('--input-prefix', dest='prefix_in', default='NONE', type=str,
                  help='Input prefix. '+
                  'Systematics maps will be input-prefix + _<syst_name>.fits. '+
                  'Mask will be input-prefix + _MaskedFraction.fits')
parser.add_option('--input-maps', dest='fname_maps', default='NONE', type=str,
                  help='Path to input maps file')
parser.add_option('--ell-bins', dest='fname_ellbins', default='NONE', type=str,
                  help='Path to ell-binning file. '+
                  'Format should be: double column (l_min,l_max). One row per bandpower.')
parser.add_option('--output-file', dest='prefix_out',default=None,type=str,
                  help='Output file name. Output will be in SACC format.')
parser.add_option('--analysis-band', dest='band', default='i', type=str,
                  help='Band considered for your analysis (g,r,i,z,y)')
parser.add_option('--mcm-output', dest='fname_mcm', default='NONE', type=str,
                  help='File containing the mode-coupling matrix. '+
                  'If NONE or non-existing, it will be computed. '+
                  'If not NONE and non-existing, new file will be created')
parser.add_option('--depth-cut', dest='depth_cut', default=24.5, type=float,
                  help='Minimum depth to consider in your footprint')
parser.add_option('--masking-threshold',dest='mask_thr',default=0.5, type=float,
                  help='Will discard all pixel with a masked fraction larger than this.')
parser.add_option('--hsc-field',dest='hsc_field',default='HSC_WIDE',type=str,
                  help="HSC field used here (just for labelling purposes)")
parser.add_option('--cont-depth',dest='cont_depth',default=False,action='store_true',
                  help='Marginalize over depth map template')
parser.add_option('--cont-dust',dest='cont_dust',default=False,action='store_true',
                  help='Marginalize over dust template')
parser.add_option('--cont-dust-bands',dest='cont_dust_bands',default=False,
                  action='store_true',help='Marginalize over dust template in all bands')
parser.add_option('--cont-stars',dest='cont_stars',default=False,action='store_true',
                  help='Marginalize over stars template')
parser.add_option('--cont-oc',dest='cont_oc',type='string',
                  action='callback',callback=opt_callback,
                  help='If you want to marginalize over observing condition maps, '+
                  'list here all observing conditions you want to include')
parser.add_option('--cont-oc-bands',dest='cont_oc_bands',default=False,action='store_true',
                  help='Marginalize over observing contition templates in all bands.')
parser.add_option('--cont-deproj-bias',dest='cont_deproj_bias',default=False,action='store_true',
                  help='Remove deprojection bias.')
parser.add_option('--covariance-option',dest='covar_opt',default='NONE',type=str,
                  help='Option to compute the covariance matrix. Options are \'NONE\''+
                  '(no covariance),  \'analytic\' (Gaussian covariance),'+
                  ' or \'gaus_sim\' (MC over Gaussian simulations)')
parser.add_option('--guess-cell',dest='guess_cell',default='NONE',type=str,
                  help='Choice of best-guess power spectra to use. Options are: \'theory\''+
                  '(input theoretical power spectra, given by --theory-prediction, '+
                  'or \'data\' (Cls from data)')
parser.add_option('--covariance-ssc',dest='compute_ssc',default=False,action='store_true',
                  help='Compute super-sample covariance part.')
parser.add_option('--theory-prediction',dest='cl_theory',default='NONE',type=str,
                  help='If computing the covariance or deprojection bias from theory, you must supply a file with a'+
                  'theory prediction for the power spectra. First column should be a'+
                  'list of ell values, subsequent columns should contain all unique cross '+
                  'correlations (11,12,...,1N,22,23,...,NN)')
parser.add_option('--covariance-coupling-file',dest='covar_coup',default='NONE',type=str,
                  help='If computing the theory covariance, pass a file name where the '+
                  'coupling coefficients will be stored. If NONE, they won\'t be saved')
parser.add_option('--syst-masking-file',dest='syst_mask_file',default='NONE',type=str,
                  help='Path to a file containing a list of systematics that one should mask. The format should be '+
                  'four columns: 1- systematic name, 2- filter (u,g,r,i,z), 3- > or <, 4- Threshold value')
parser.add_option('--noise-bias',dest='noise_bias',default='analytic',type=str,
                  help='Option to compute the noise bias. Options are: \'analytic\''+
                  'or \'simulated\'.')
parser.add_option('--nrealiz',dest='nrealiz',type=int,default=1000,
                  help='If the noise bias is computed using simulations, this parameters determines the number of'+
                  'simulations to use.')
parser.add_option('--bo-mask-type', dest='mask_type', default='sirius', type=str,
                  help='Bright object mask (arcturus or sirius)')

####
# Read options
(o, args) = parser.parse_args()

prefix_in_use=o.prefix_in
if o.mask_type=='arcturus' :
  prefix_in_use+='_marct'

def read_map_bands(fname,read_bands,bandname) :
  if read_bands :
    i_map=-1
  else :
    i_map=['g','r','i','z','y'].index(bandname)
  fskb,temp=fm.read_flat_map(fname,i_map=i_map)
  fm.compare_infos(fsk,fskb)
  if i_map!=-1 :
    temp=[temp]
    
  return temp


print("Reading mask")
#Create depth-based mask
fsk,mp_depth=fm.read_flat_map(prefix_in_use+"_10s_depth_mean_fluxerr.fits",2)
mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=o.depth_cut]=1

#Read masked fraction
fskb,mskfrac=fm.read_flat_map(prefix_in_use+'_MaskedFraction.fits',i_map=0)
fm.compare_infos(fsk,fskb)

#Create BO-based mask
msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>o.mask_thr]=1
msk_t=msk_bo*msk_depth

#Mask systematics
do_mask_syst=not (o.syst_mask_file=='NONE')
msk_syst=msk_t.copy()
if do_mask_syst :
  #Read systematics cut data
  data_syst=np.genfromtxt(o.syst_mask_file,dtype=[('name','|U32'),('band','|U4'),('gl','|U4'),('thr','<f8')])
  for d in data_syst :
    #Read systematic
    if d['name'].startswith('oc_') :
      sysmap=read_map_bands(o.prefix_in+'_'+d['name']+'.fits',False,d['band'])[0]
    elif d['name']=='dust' :
      sysmap=read_map_bands(o.prefix_in+'_syst_dust.fits',False,d['band'])[0]
    else :
      raise KeyError("Unknown systematic name "+d['name'])
    
    #Divide by mean
    sysmean=np.sum(msk_t*mskfrac*sysmap)/np.sum(msk_t*mskfrac)
    sysmap/=sysmean

    #Apply threshold
    msk_sys_this=msk_t.copy(); fsky_pre=np.sum(msk_syst)
    if d['gl']=='<' :
      msk_sys_this[sysmap<d['thr']]=0
    else :
      msk_sys_this[sysmap>d['thr']]=0
    msk_syst*=msk_sys_this
    fsky_post=np.sum(msk_syst)
    print(' '+d['name']+d['gl']+'%.3lf'%(d['thr'])+
          ' removes ~%.2lf per-cent of the available sky'%((1-fsky_post/fsky_pre)*100))
  print(' All systematics remove %.2lf per-cent of the sky'%((1-np.sum(msk_syst)/np.sum(msk_t))*100))
msk_t*=msk_syst

#Area
area_patch=np.sum(msk_t*mskfrac)*np.radians(fsk.dx)*np.radians(fsk.dy)

#Read contaminants
print("Reading contaminant templates")
temps=[]
# 1- Depth
if o.cont_depth :
  temps.append(mp_depth)
# 2- Dust
if o.cont_dust :
  temp=read_map_bands(o.prefix_in+'_syst_dust.fits',o.cont_dust_bands,o.band)
  for t in temp :
    temps.append(t)
# 3- Stars
if o.cont_stars :
  fskb,temp=fm.read_flat_map(o.prefix_in+'_syst_nstar_'+o.band+'%.2lf.fits'%o.depth_cut,
                             i_map=0)
  fm.compare_infos(fsk,fskb)
  temps.append(temp)
# 4- Observing conditions
if o.cont_oc is not None :
  oc_all=['airmass','ccdtemp','ellipt','exptime','nvisit','seeing','sigma_sky','skylevel']
  #Include only OCs we know about
  oc_list=[]
  for c in o.cont_oc :
      if c in oc_all :
        oc_list.append(c)
  for c in oc_list :
    temp=read_map_bands(o.prefix_in+'_oc_'+c+'.fits',o.cont_oc_bands,o.band)
    for t in temp :
      temps.append(t)
if temps==[] :
  temps=None
else :
  print(" - Will marginalize over a total of %d contaminant templates"%(len(temps)))
  #Remove mean
  for i_t,t in enumerate(temps) :
    temps[i_t]-=np.sum(msk_t*mskfrac*t)/np.sum(msk_t*mskfrac)

#Set binning scheme
lini,lend=np.loadtxt(o.fname_ellbins,unpack=True)
bpws=nmt.NmtBinFlat(lini,lend)
ell_eff=bpws.get_effective_ells()

#Generate tracers
print("Generating tracers")
class Tracer(object) :
  def __init__(self,hdu_list,i_bin,fsk,mask_binary,masked_fraction,contaminants=None) :
    
    #Read numbers map
    self.fsk,nmap=fm.read_flat_map(None,hdu=hdu_list[2*i_bin])

    #Read N(z)
    self.nz_data=hdu_list[2*i_bin+1].data.copy()

    #Make sure other maps are compatible
    fm.compare_infos(self.fsk,fsk)
    if not self.fsk.is_map_compatible(mask_binary) :
      raise ValueError("Mask size is incompatible")
    if not self.fsk.is_map_compatible(masked_fraction) :
      raise ValueError("Mask size is incompatible")
    if contaminants is not None :
      for ic,c in enumerate(contaminants) :
        if not self.fsk.is_map_compatible(c) :
          raise ValueError("%d-th contaminant template is incompatible"%ic)
          
    #Translate into delta map
    self.masked_fraction=masked_fraction
    self.weight=masked_fraction*mask_binary
    goodpix=np.where(mask_binary>0.1)[0]
    self.goodpix=goodpix
    self.mask_binary=mask_binary
    self.Ngal = np.sum(nmap*mask_binary)
    ndens=np.sum(nmap*mask_binary)/np.sum(self.weight)
    self.ndens_perad=ndens/(np.radians(self.fsk.dx)*np.radians(self.fsk.dy))
    self.delta=np.zeros_like(self.weight)
    self.delta[goodpix]=nmap[goodpix]/(ndens*masked_fraction[goodpix])-1

    #Reshape contaminants
    conts=None
    if contaminants is not None :
      conts=[[c.reshape([self.fsk.ny,self.fsk.nx])] for c in contaminants]

    #Form NaMaster field
    self.field=nmt.NmtFieldFlat(np.radians(self.fsk.lx),np.radians(self.fsk.ly),
                                self.weight.reshape([self.fsk.ny,self.fsk.nx]),
                                [self.delta.reshape([self.fsk.ny,self.fsk.nx])],
                                templates=conts)

hdul=fits.open(o.fname_maps)
if len(hdul)%2!=0 :
  raise ValueError("Input file should have two HDUs per map")
nbins=len(hdul)//2
tracers=[Tracer(hdul,i,fsk,msk_t,mskfrac,contaminants=temps) for i in range(nbins)]
hdul.close()

#Compute MCM
wsp=nmt.NmtWorkspaceFlat()
if not os.path.isfile(o.fname_mcm) :
  print("Computing mode-coupling matrix")
  wsp.compute_coupling_matrix(tracers[0].field,tracers[0].field,bpws)
  if o.fname_mcm!='NONE' :
    wsp.write_to(o.fname_mcm)
else :
  print("Reading mode-coupling matrix from file")
  wsp.read_from(o.fname_mcm)

#Compute window functions
lmax=int(180.*np.sqrt(1./fsk.dx**2+1./fsk.dy**2))
nbands=wsp.wsp.bin.n_bands
l_arr=np.arange(lmax+1)
if (o.fname_mcm=='NONE') or (not os.path.isfile(o.fname_mcm+".windows.npz")) :
  print("Computing window functions")
  windows=np.zeros([nbands,lmax+1])
  t_hat=np.zeros(lmax+1);
  for il,l in enumerate(l_arr) :
    t_hat[il]=1.;
    windows[:,il]=wsp.decouple_cell(wsp.couple_cell(l_arr,[t_hat]))
    t_hat[il]=0.;
  if o.fname_mcm!='NONE' :
    np.savez(o.fname_mcm+".windows",windows=windows)
else :
  print("Reading window functions")
  windows=np.load(o.fname_mcm+".windows.npz")['windows']
    

#Compute all cross-correlations
print("Computing all cross-power specra")
cls_all=[]
windows_sacc=[]
ordering=np.zeros([nbins,nbins],dtype=int)
i_x=0
for i in range(nbins) :
  for j in range(i,nbins) :
    cl_coupled=nmt.compute_coupled_cell_flat(tracers[i].field,tracers[j].field,bpws)
    cls_all.append(wsp.decouple_cell(cl_coupled)[0])
    for b in range(nbands) :
      windows_sacc.append(sacc.Window(l_arr,windows[b]))
    ordering[i,j]=i_x
    if j!=i :
      ordering[j,i]=i_x
    i_x+=1
cls_all=np.array(cls_all)
n_cross=len(cls_all)
n_ell=len(ell_eff)

#Get guess power spectra
if o.guess_cell=='data' :
  #Interpolate measured power spectra
  lth=np.arange(2,lmax+1)+0.
  clth=np.zeros([n_cross,len(lth)])
  for i in range(n_cross) :
    clf=interp1d(ell_eff,cls_all[i],bounds_error=False,fill_value=0,kind='linear')
    clth[i,:]=clf(lth)
    clth[i,lth<=ell_eff[0]]=cls_all[i,0]
    clth[i,lth>=ell_eff[-1]]=cls_all[i,-1]
elif o.guess_cell=='theory' :
  #Read theory power spectra
  data=np.loadtxt(o.cl_theory,unpack=True)
  lth=data[0]
  clth=data[1:]
else :
  raise ValueError("Must provide a valid guess C_ell\n")
if len(clth)!=n_cross :
  raise ValueError("Theory power spectra have a wrong shape")
    
#Compute deprojection bias
cl_deproj_all=[]
if o.cont_deproj_bias :
  print("Computing deprojection bias")
  cls_all=[]
i_x=0
for i in range(nbins) :
  for j in range(i,nbins) :
    if o.cont_deproj_bias :
      cl_coupled=nmt.compute_coupled_cell_flat(tracers[i].field,tracers[j].field,bpws)
      cl_deproj_bias=nmt.deprojection_bias_flat(tracers[i].field,tracers[j].field,bpws,lth,[clth[i_x]])
      cls_all.append(wsp.decouple_cell(cl_coupled,cl_bias=cl_deproj_bias)[0])
    else :
      cl_deproj_bias=None
    cl_deproj_all.append(cl_deproj_bias)
    i_x+=1
if o.cont_deproj_bias :
  cls_all=np.array(cls_all)

#Compute covariance matrix
if o.covar_opt=='NONE' :
  covar=None
elif o.covar_opt=='analytic' :
  print("Computing analytic Gaussian covariance")
  #Initialize covariance
  covar=np.zeros([n_cross*n_ell,n_cross*n_ell])

  #Compute coupling coefficients
  cwsp=nmt.NmtCovarianceWorkspaceFlat();
  if not os.path.isfile(o.covar_coup) :
    cwsp.compute_coupling_coefficients(wsp,wsp)
    if o.covar_coup!='NONE' :
      cwsp.write_to(o.covar_coup)
  else :
    cwsp.read_from(o.covar_coup)

  ix_1=0
  for i1 in range(nbins) :
    for j1 in range(i1,nbins) :
      ix_2=0
      for i2 in range(nbins) :
        for j2 in range(i2,nbins) :
          ca1b1=clth[ordering[i1,i2]]
          ca1b2=clth[ordering[i1,j2]]
          ca2b1=clth[ordering[j1,i2]]
          ca2b2=clth[ordering[j1,j2]]
          cov_here=nmt.gaussian_covariance_flat(cwsp,lth,ca1b1,ca1b2,ca2b1,ca2b2)
          covar[ix_1*n_ell:(ix_1+1)*n_ell,:][:,ix_2*n_ell:(ix_2+1)*n_ell]=cov_here
          ix_2+=1
      ix_1+=1

elif o.covar_opt=='gaus_sim' :
  #Read theory power spectra
  nsims=10*n_cross*n_ell #Use 10 times as many simulations as there are data points
  print("Computing covariance from %d Gaussian simulations"%nsims)
  msk_binary=msk_t.reshape([fsk.ny,fsk.nx])
  weights=(msk_t*mskfrac).reshape([fsk.ny,fsk.nx])
  if temps is not None :
    conts=[[t.reshape([fsk.ny,fsk.nx])] for t in temps]
  else :
    conts=None
  
  cells_sims=[]
  for isim in range(nsims) :
    if isim%100==0 :
      print(" %d-th isim"%isim)
    #Generate random maps
    mps=nmt.synfast_flat(fsk.nx,fsk.ny,np.radians(fsk.lx),np.radians(fsk.ly),clth,np.zeros(nbins),seed=1000+isim)
    #Nmt fields
    flds=[nmt.NmtFieldFlat(np.radians(fsk.lx),np.radians(fsk.ly),weights,[m],templates=conts) for m in mps]
    #Compute power spectra (possibly with deprojection)
    i_x=0
    cells_this=[]
    for i in range(nbins) :
      for j in range(i,nbins) :
        cells_this.append(wsp.decouple_cell(nmt.compute_coupled_cell_flat(flds[i],flds[j],bpws),
                                            cl_bias=cl_deproj_all[i_x])[0])
        i_x+=1
    cells_sims.append(np.array(cells_this).flatten())
  cells_sims=np.array(cells_sims)
  #Save simulations for further analysis
  np.savez(o.prefix_out+".cls",cl_sims=cells_sims)

  #Compute covariance
  covar=np.cov(cells_sims.T)
else :
  print("Unknown covariance option "+o.covar_opt+" no covariance computed")
  covar=None

if (covar is not None) and o.compute_ssc :
  #Compute number density and uncertainty on it from cosmic variance
  import pyccl as ccl
  from scipy.special import jv
  #Cosmology
  cosmo=ccl.Cosmology(Omega_c=0.27,Omega_b=0.049,h=0.67,sigma8=0.83,w0=-1.,wa=0.,n_s=0.96)

  #Sky fraction
  f_sky=np.sum(msk_t*mskfrac)*fsk.dx*fsk.dy/(4*np.pi*(180./np.pi)**2)

  #Tracers
  cclt=[]
  z_b=np.array([0.0,0.5,1.0,2.0,4.0]);
  b_b=np.array([0.82,1.10,1.44,1.66,2.61])
  b_bf=interp1d(z_b,b_b)
  ng_data=np.zeros([len(tracers),4]);
  for i_t,t in enumerate(tracers) :
    zarr=(t.nz_data['z_i']+t.nz_data['z_f'])*0.5
    narr=t.nz_data['n_z']
    barr=b_bf(zarr)
    cclt.append(ccl.NumberCountsTracer(cosmo,has_rsd=False,dndz=(zarr,narr),
                                       bias=(zarr,barr)))

    lmax_ng=3000
    larr=np.arange(lmax_ng+1)
    cell=ccl.angular_cl(cosmo,cclt[i_t],cclt[i_t],larr) # A = pi*th^2
    theta_s=np.sqrt(area_patch/np.pi)
    well=np.ones(len(larr)); well[1:]=2*jv(1,larr[1:]*theta_s)/(larr[1:]*theta_s); well=well**2
    ngals=t.ndens_perad*(np.radians(t.fsk.dx)*np.radians(t.fsk.dy))*np.sum(t.weight)
    sigma_c=ngals*np.sqrt(np.sum(larr*cell*well/(2*np.pi)))
    sigma_p=np.sqrt(ngals)
    ng_data[i_t,0]=ngals/area_patch
    ng_data[i_t,1]=sigma_p/area_patch
    ng_data[i_t,2]=sigma_c/area_patch
    ng_data[i_t,3]=np.sqrt(sigma_p**2+sigma_c**2)/area_patch

  print(ng_data)
  np.savetxt(o.prefix_out+".ngals",ng_data,header='[1]-Ngals [2]-Sigma_poisson [3]-Sigma_CV [4]-Sigma_T')

  #SSC init
  def get_response_func() :
    zarr=np.array([4,3,2,1,0])
    resp2d=[]
    for iz,z in enumerate(zarr) :
      kresph,_,_,_,resp1d=np.loadtxt("ssc_responses/Response_z%d.txt"%(int(z)),unpack=True)
      resp2d.append(resp1d)
    kresp=kresph*0.67
    aresp=1./(1.+zarr)
    resp2d=np.array(resp2d)
    return ccl.Pk2D(a_arr=aresp,lk_arr=np.log(kresp),pk_arr=resp2d,is_logp=False)
  respf=get_response_func()
  ssc_wsp=ccl.SSCWorkspace(cosmo,f_sky,cclt[0],cclt[0],cclt[0],cclt[0],ell_eff,respf)

  #SSC compute
  covar_ssc=np.zeros([n_cross*n_ell,n_cross*n_ell])
  iv=0
  for i1 in range(nbins) :
    for i2  in range(i1,nbins) :
      jv=0
      for j1 in range(nbins) :
        for j2  in range(j1,nbins) :
          mat=ccl.angular_cl_ssc_from_workspace(ssc_wsp,cosmo,
                                                cclt[i1],cclt[i2],
                                                cclt[j1],cclt[j2])
          covar_ssc[iv*n_ell:(iv+1)*n_ell,jv*n_ell:(jv+1)*n_ell]=mat
          jv+=1
      iv+=1

  covar+=covar_ssc

#Compute noise bias
if o.noise_bias == 'analytic':
    nls_all=np.zeros_like(cls_all)
    i_x=0
    for i in range(nbins) :
      for j in range(i,nbins) :
        if i==j : #Add shot noise in auto-correlation
          t=tracers[i]
          corrfac=np.sum(t.weight)/(t.fsk.nx*t.fsk.ny)
          nl=np.ones_like(ell_eff)*corrfac/t.ndens_perad
          nls_all[i_x]=wsp.decouple_cell([nl])[0]
        i_x+=1

elif o.noise_bias == 'simulated':

    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.info('Computing shape noise power spectrum from Poisson simulations.')
    logger.info('Running nrealiz = {}.'.format(o.nrealiz))

    def randomize_deltag_map(tracer):
        """
        Creates a randomised version of the input map map by assigning the
        galaxies in the surevy to random pixels in the map. Basically it rotates each
        galaxy by a random angle but not rotating it out of the survey footprint.
        :param map: masked galaxy overdensity map which needs to randomised
        :param Ngal: number of galaxies used to create the map
        :return randomised_map: a randomised version of the masked input map
        """

        logger.info('Randomizing galaxy map.')

        mask = tracer.weight.reshape([tracer.fsk.ny, tracer.fsk.nx])
        Ngal = int(tracer.Ngal)

        np.random.seed(seed=None)
        maskpixy, maskpixx = np.where(mask != 0.)

        galpix_mask = np.random.choice(np.arange(maskpixx.shape[0]), size=Ngal, \
                                       p=mask[mask != 0.]/np.sum(mask[mask != 0.]))

        galpixx = maskpixx[galpix_mask]
        galpixy = maskpixy[galpix_mask]

        maskshape = mask.shape
        ny, nx = maskshape
        ipix = galpixx + nx*galpixy

        randomized_nmap = np.bincount(ipix, minlength=nx*ny)

        randomized_deltamap = np.zeros_like(randomized_nmap, dtype='float')
        ndens = np.sum(randomized_nmap*tracer.mask_binary)/np.sum(tracer.weight)
        randomized_deltamap[tracer.goodpix] = randomized_nmap[tracer.goodpix]/(ndens*tracer.masked_fraction[tracer.goodpix])-1
        randomized_deltamap = randomized_deltamap.reshape(maskshape)

        return randomized_deltamap

    nls_all = np.zeros_like(cls_all)
    wsps = [[None for i in range(nbins)] for ii in range(nbins)]
    i_x = 0
    for i in range(nbins) :
      for j in range(i,nbins) :
        if i == j: #Add shot noise in auto-correlation
            tracer = tracers[i]
            mask = tracer.weight.reshape([tracer.fsk.ny, tracer.fsk.nx])
            ncl_uncoupled = np.zeros((o.nrealiz, ell_eff.shape[0]))
            for ii in range(o.nrealiz):
                randomized_map = randomize_deltag_map(tracer)
                f0 = nmt.NmtFieldFlat(np.radians(fsk.lx),np.radians(fsk.ly), mask, [randomized_map], purify_b=False)
                # if wsps[i][j] == None:
                #     logger.info('Workspace element for i, j = {}, {} not set.'.format(i, j))
                #     logger.info('Computing workspace element.')
                #     wsp = nmt.NmtWorkspaceFlat()
                #     wsp.compute_coupling_matrix(f0, f0, bpws)
                #     wsps[i][j] = wsp
                # else:
                #     logger.info('Workspace element already set for i, j = {}, {}.'.format(i, j))

                # Compute pseudo-Cls
                ncl_coupled = nmt.compute_coupled_cell_flat(f0, f0, bpws)
                # Uncoupling pseudo-Cls
                ncl_uncoupled[ii, :] = wsp.decouple_cell(ncl_coupled)
            ncl_uncoupled_mean = np.mean(ncl_uncoupled, axis=0)
            nls_all[i_x] = ncl_uncoupled_mean
        i_x+=1

#Save to SACC format
print("Saving to SACC")
#Tracers
sacc_tracers=[]
for i_t,t in enumerate(tracers) :
  z=(t.nz_data['z_i']+t.nz_data['z_f'])*0.5
  nz=t.nz_data['n_z']
  T=sacc.Tracer('bin_%d'%i_t,'point',z,nz,exp_sample=o.hsc_field)
  T.addColumns({'ndens':t.ndens_perad*np.ones_like(nz)})
  sacc_tracers.append(T)
#Binning and mean
type,ell,dell,t1,q1,t2,q2=[],[],[],[],[],[],[]
for t1i in range(nbins) :
  for t2i in range(t1i,nbins) :
    for i_l,l in enumerate(ell_eff) :
      type.append('F') #Fourier-space
      ell.append(l)
      dell.append(lend[i_l]-lini[i_l])
      t1.append(t1i)
      q1.append('P')
      t2.append(t2i)
      q2.append('P')
sacc_binning=sacc.Binning(type,ell,t1,q1,t2,q2,deltaLS=dell,windows=windows_sacc)
sacc_mean=sacc.MeanVec(cls_all.flatten())
if covar is None :
  sacc_precision=None
else :
  sacc_precision=sacc.Precision(covar,"dense",is_covariance=True, binning=sacc_binning)
sacc_meta={'Field':o.hsc_field,'Area_rad':area_patch}
s=sacc.SACC(sacc_tracers,sacc_binning,sacc_mean,precision=sacc_precision,meta=sacc_meta)
s.printInfo()
s.saveToHDF(o.prefix_out+'.sacc')
#Save noise
sacc_binning_noise=sacc.Binning(type,ell,t1,q1,t2,q2,deltaLS=dell)
sacc_mean_noise=sacc.MeanVec(nls_all.flatten())
s=sacc.SACC(sacc_tracers,sacc_binning_noise,sacc_mean_noise,precision=None,meta=sacc_meta)
s.printInfo()
s.saveToHDF(o.prefix_out+'_noise.sacc')
