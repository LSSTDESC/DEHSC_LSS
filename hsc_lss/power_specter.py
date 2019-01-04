from ceci import PipelineStage
from .types import FitsFile,ASCIIFile,BinaryFile,NpzFile,SACCFile
import numpy as np
from .flatmaps import read_flat_map,compare_infos
from astropy.io import fits
import pymaster as nmt
from .tracer import Tracer
import os
import sacc
from scipy.interpolate import interp1d

class PowerSpecter(PipelineStage) :
    name="PowerSpecter"
    inputs=[('masked_fraction',FitsFile),('ngal_maps',FitsFile),
            ('dust_map',FitsFile),('star_map',FitsFile),('depth_map',FitsFile),
            ('ccdtemp_maps',FitsFile),('airmass_maps',FitsFile),('exptime_maps',FitsFile),
            ('skylevel_maps',FitsFile),('sigma_sky_maps',FitsFile),('seeing_maps',FitsFile),
            ('ellipt_maps',FitsFile),('nvisit_maps',FitsFile),('cosmos_weights',FitsFile)]
    outputs=[('mcm',BinaryFile),('cov_mcm',BinaryFile),('gaucov_sims',NpzFile),
             ('windows_l',NpzFile),('noi_bias',SACCFile),('dpj_bias',SACCFile),
             ('power_spectra_wdpj',SACCFile),('power_spectra_wodpj',SACCFile)]
    config_options={'ell_bpws':[100.0,200.0,300.0,
                                400.0,600.0,800.0,
                                1000.0,1400.0,1800.0,
                                2200.0,3000.0,3800.0,
                                4600.0,6200.0,7800.0,
                                9400.0,12600.0,12600.0,15800.0],
                    'z_bias_nodes':[0.00,0.50,1.00,2.00,4.00],
                    'b_bias_nodes':[0.82,1.10,1.44,1.66,2.61],
                    'depth_cut':24.5,'band':'i','mask_thr':0.5,'guess_spectrum':'NONE',
                    'gaus_covar_type':'analytic','oc_all_bands':True,'add_ssc':False,
                    'mask_systematics':False,'noise_bias_type':'analytic'}

    def read_map_bands(self,fname,read_bands,bandname) :
        if read_bands :
            i_map=-1
        else :
            i_map=['g','r','i','z','y'].index(bandname)
        fskb,temp=read_flat_map(fname,i_map=i_map)
        compare_infos(self.fsk,fskb)
        if i_map!=-1 :
            temp=[temp]
    
        return temp

    def get_sacc_windows(self,wsp) :
        #Compute window functions
        nbands=wsp.wsp.bin.n_bands
        l_arr=np.arange(self.lmax+1)
        if not os.path.isfile(self.get_output('windows_l')) :
            print("Computing window functions")
            windows=np.zeros([nbands,self.lmax+1])
            t_hat=np.zeros(self.lmax+1);
            for il,l in enumerate(l_arr) :
                t_hat[il]=1.;
                windows[:,il]=wsp.decouple_cell(wsp.couple_cell(l_arr,[t_hat]))
                t_hat[il]=0.;
            np.savez(self.get_output('windows_l')[:-4],windows=windows)
        else :
            print("Reading window functions")
            windows=np.load(self.get_output('windows_l'))['windows']

        windows_sacc=[]
        #i_x=0
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                for b in range(nbands) :
                    windows_sacc.append(sacc.Window(l_arr,windows[b]))

        return windows_sacc

    def get_noise(self,tracers,wsp,bpws,nsims=1000) :
        if self.config['noise_bias_type']=='analytic' :
            return self.get_noise_analytic(tracers,wsp)
        elif self.config['noise_bias_type']=='pois_sim' :
            return self.get_noise_simulated(tracers,wsp,bpws,nsims)

    def get_noise_analytic(self,tracers,wsp) :
        nls_all=np.zeros([self.ncross,self.nell])
        i_x=0
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                if i==j : #Add shot noise in auto-correlation
                    t=tracers[i]
                    corrfac=np.sum(t.weight)/(t.fsk.nx*t.fsk.ny)
                    nl=np.ones(self.nell)*corrfac/t.ndens_perad
                    nls_all[i_x]=wsp.decouple_cell([nl])[0]
                i_x+=1
        return nls_all
        
    def get_noise_simulated(self,tracers,wsp,bpws,nsims) :
        def randomize_delgag_map(tracer,seed) :
            """
            Creates a randomised version of the input map map by assigning the
            galaxies in the surevy to random pixels in the map. Basically it rotates each
            galaxy by a random angle but not rotating it out of the survey footprint.
            :param map: masked galaxy overdensity map which needs to randomised
            :param Ngal: number of galaxies used to create the map
            :return randomised_map: a randomised version of the masked input map
            """
            
            mask = tracer.weight.reshape([tracer.fsk.ny, tracer.fsk.nx])
            Ngal = int(tracer.Ngal)

            np.random.seed(seed=seed)
            maskpixy,maskpixx=np.where(mask!=0.)
            galpix_mask=np.random.choice(np.arange(maskpixx.shape[0]),size=Ngal,
                                         p=mask[mask != 0.]/np.sum(mask[mask != 0.]))
            galpixx=maskpixx[galpix_mask]
            galpixy=maskpixy[galpix_mask]

            maskshape=mask.shape
            ny,nx=maskshape
            ipix=galpixx+nx*galpixy

            randomized_nmap=np.bincount(ipix,minlength=nx*ny)

            randomized_deltamap=np.zeros_like(randomized_nmap,dtype='float')
            ndens=np.sum(randomized_nmap*tracer.mask_binary)/np.sum(tracer.weight)
            randomized_deltamap[tracer.goodpix]=randomized_nmap[tracer.goodpix]/(ndens*tracer.masked_fraction[tracer.goodpix])-1
            randomized_deltamap=randomized_deltamap.reshape(maskshape)

            return randomized_deltamap

        nls_all=np.zeros([self.ncross,self.nell])
        wsps=[[None for i in range(self.nbins)] for ii in range[nbins]]
        i_x=0
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                if i==j: #Add shot noise in the auto-correlation
                    tracer=tracers[i]
                    mask=tracer.weight.reshape([tracer.fsk.ny,tracer.fsk.nx])
                    ncl_uncoupled=np.zeros((nsims,self.nell))
                    for ii in range(nsims) :
                        randomized_map=randomize_deltag_map(tracer)
                        f0=nmt.NmtFieldFlat(np.radians(self.fsk.lx),np.radians(self.fsk.ly),mask,
                                            [randomized_map])
                        ncl_uncoupled[ii,:]=wsp.decouple_cell(nmt.coupled_cell_flat(f0,f0,bpws))
                    nls_all[i_x]=np.mean(ncl_uncoupled,axis=0)
                i_x+=1

        return nls_all

    def get_covar_ssc(self,tracers,ell_eff) :
        #Compute number density and uncertainty on it from cosmic variance
        import pyccl as ccl
        from scipy.special import jv
        #Cosmology
        cosmo=ccl.Cosmology(Omega_c=0.27,Omega_b=0.049,h=0.67,sigma8=0.83,w0=-1.,wa=0.,n_s=0.96)

        #Sky fraction
        f_sky=np.sum(self.msk_bi*self.mskfrac)*self.area_pix/(4*np.pi)

        #Tracers
        cclt=[]
        z_b=np.array(self.config['z_bias_nodes'])
        b_b=np.array(self.config['b_bias_nodes'])
        b_bf=interp1d(z_b,b_b)
        ng_data=np.zeros([len(tracers),4]);
        for i_t,t in enumerate(tracers) :
            zarr=(t.nz_data['z_i']+t.nz_data['z_f'])*0.5
            narr=t.nz_data['nz_cosmos']
            barr=b_bf(zarr)
            cclt.append(ccl.NumberCountsTracer(cosmo,has_rsd=False,dndz=(zarr,narr),
                                               bias=(zarr,barr)))

            larr=np.arange(3001)
            cell=ccl.angular_cl(cosmo,cclt[i_t],cclt[i_t],larr) # A = pi*th^2
            theta_s=np.sqrt(self.area_patch/np.pi)
            well=np.ones(len(larr)); well[1:]=2*jv(1,larr[1:]*theta_s)/(larr[1:]*theta_s); well=well**2
            ngals=t.ndens_perad*(np.radians(t.fsk.dx)*np.radians(t.fsk.dy))*np.sum(t.weight)
            sigma_c=ngals*np.sqrt(np.sum(larr*cell*well/(2*np.pi)))
            sigma_p=np.sqrt(ngals)
            ng_data[i_t,0]=ngals/self.area_patch
            ng_data[i_t,1]=sigma_p/self.area_patch
            ng_data[i_t,2]=sigma_c/self.area_patch
            ng_data[i_t,3]=np.sqrt(sigma_p**2+sigma_c**2)/self.area_patch

        #SSC init
        def get_response_func() :
            zarr=np.array([4.,3.,2.,1.,0.])
            resp2d=[]
            for iz,z in enumerate(zarr) :
                kresph,_,_,_,resp1d=np.loadtxt(self.config['ssc_response_prefix']+"_z%d.txt"%(int(z)),
                                               unpack=True)
                resp2d.append(resp1d)
            kresp=kresph*0.67
            aresp=1./(1.+zarr)
            resp2d=np.array(resp2d)
            return ccl.Pk2D(a_arr=aresp,lk_arr=np.log(kresp),pk_arr=resp2d,is_logp=False)
        respf=get_response_func()
        ssc_wsp=ccl.SSCWorkspace(cosmo,f_sky,cclt[0],cclt[0],cclt[0],cclt[0],ell_eff,respf)

        #SSC compute
        covar_ssc=np.zeros([self.ncross*self.nell,self.ncross*self.nell])
        iv=0
        for i1 in range(self.nbins) :
            for i2  in range(i1,self.nbins) :
                jv=0
                for j1 in range(self.nbins) :
                    for j2  in range(j1,self.nbins) :
                        mat=ccl.angular_cl_ssc_from_workspace(ssc_wsp,cosmo,
                                                              cclt[i1],cclt[i2],
                                                              cclt[j1],cclt[j2])
                        covar_ssc[iv*self.nell:(iv+1)*self.nell,jv*self.nell:(jv+1)*self.nell]=mat
                        jv+=1
                iv+=1
        
        return ng_data,covar_ssc

    def get_dpj_bias(self,trc,lth,clth,cl_coupled,wsp,bpws) :
        #Compute deprojection bias
        if os.path.isfile(self.get_output('dpj_bias')) :
            print("Reading deprojection bias")
            s=sacc.SACC.loadFromHDF(self.get_output('dpj_bias'))
            cls_deproj_all=s.mean.vector.reshape([self.ncross,self.nell])
        else :
            print("Computing deprojection bias")
            cls_deproj_all=[]
            i_x=0
            for i in range(self.nbins) :
                for j in range(i,self.nbins) :
                    print(i,j)
                    cl_deproj_bias=nmt.deprojection_bias_flat(trc[i].field,trc[j].field,bpws,
                                                              lth,[clth[i_x]])[0]
                    cls_deproj_all.append(cl_deproj_bias)
                    i_x+=1
            cls_deproj_all=np.array(cls_deproj_all)

        #Remove deprojection bias
        cls_all=[]
        i_x=0
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                cls_all.append(wsp.decouple_cell([cl_coupled[i_x]],
                                                 cl_bias=[cls_deproj_all[i_x]])[0])
                i_x+=1
        return np.array(cls_all),cls_deproj_all

    def get_cl_guess(self,ld,cld) :
        if self.config['guess_spectrum']=='NONE' :
            print("Interpolating data power spectra")
            lth=np.arange(2,self.lmax+1)
            clth=np.zeros([self.ncross,len(lth)])
            for i in range(self.ncross) :
                clf=interp1d(ld,cld[i],bounds_error=False,fill_value=0,kind='linear')
                clth[i,:]=clf(lth)
                clth[i,lth<=ld[0]]=cld[i,0]
                clth[i,lth>=ld[-1]]=cld[i,-1]
        else :
            print("Reading theory power spectra")
            data=np.loadtxt(self.config['guess_spectrum'],unpack=True)
            lth=data[0]
            clth=data[1:]
            if len(clth)!=self.ncross :
                raise ValueError("Theory power spectra have a wrong shape")
        return lth,clth

    def get_power_spectra(self,trc,wsp,bpws) :
        cls_all=[]
        cls_coupled=[]
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                cl_coupled=nmt.compute_coupled_cell_flat(trc[i].field,trc[j].field,bpws)
                cls_all.append(wsp.decouple_cell(cl_coupled)[0])
                cls_coupled.append(cl_coupled[0])
        return np.array(cls_all),np.array(cls_coupled)

    def get_covar(self,lth,clth,bpws,wsp,temps,cl_dpj_all) :
        if self.config['gaus_covar_type']=='analytic' :
            print("Computing analytical Gaussian covariance")
            cov=self.get_covar_analytic(lth,clth,bpws,wsp)
        elif self.config['gaus_covar_type']=='gaus_sim' :
            print("Computing simulated Gaussian covariance")
            cov=self.get_covar_gaussim(lth,clth,bpws,wsp,temps,cl_dpj_all)

        return cov
            
    def get_mcm(self,tracers,bpws) :
        wsp=nmt.NmtWorkspaceFlat()
        if not os.path.isfile(self.get_output('mcm')) :
            print("Computing MCM")
            wsp.compute_coupling_matrix(tracers[0].field,tracers[0].field,bpws)
            wsp.write_to(self.get_output('mcm'))
        else :
            print("Reading MCM")
            wsp.read_from(self.get_output('mcm'))
        
        return wsp

    def get_covar_mcm(self,wsp) :
        cwsp=nmt.NmtCovarianceWorkspaceFlat()
        if not os.path.isfile(self.get_output('cov_mcm')) :
            print("Computing covariance MCM")
            cwsp.compute_coupling_coefficients(wsp,wsp)
            cwsp.write_to(self.get_output('cov_mcm'))
        else :
            print("Reading covariance MCM")
            cwsp.read_from(self.get_output('cov_mcm'))
        
        return cwsp

    def get_covar_gaussim(self,lth,clth,bpws,wsp,temps,cl_dpj_all) :
        #Create a dummy file for the covariance MCM
        f=open(self.get_output('cov_mcm'),"w")
        f.close()

        #Setup
        nsims=10*self.ncross*self.nell
        print("Computing covariance from %d Gaussian simulations"%nsims)
        msk_binary=self.msk_bi.reshape([self.fsk.ny,self.fsk.nx])
        weights=(self.msk_bi*self.mskfrac).reshape([self.fsk.ny,self.fsk.nx])
        if temps is not None :
            conts=[[t.reshape([self.fsk.ny,self.fsk.nx])] for t in temps]
            cl_dpj=[None for i in range(self.ncross)]
        else :
            conts=None
            cl_dpj=cl_dpj_all

        #Iterate
        cells_sims=[]
        for isim in range(nsims) :
            if isim%100==0 :
                print(" %d-th sim"%isim)
            #Generate random maps
            mps=nmt.synfast_flat(self.fsk.nx,self.fsk.ny,
                                 np.radians(self.fsk.lx),np.radians(self.fsk.ly),
                                 clth,np.zeros(self.nbins),seed=1000+isim)
            #Nmt fields
            flds=[nmt.NmtFieldFlat(np.radians(self.fsk.lx),np.radians(self.fsk.ly),weights,
                                   [m],templates=conts) for m in mps]
            #Compute power spectra (possibly with deprojection)
            i_x=0
            cells_this=[]
            for i in range(self.nbins) :
                for j in range(i,self.nbins) :
                    print(cl_dpj)
                    cl=nmt.compute_coupled_cell_flat(flds[i],flds[j],bpws)
                    cells_this.append(wsp.decouple_cell(cl,cl_bias=cl_dpj[i_x])[0])
                    i_x+=1
            cells_sims.append(np.array(cells_this).flatten())
        cells_sims=np.array(cells_sims)
        #Save simulations for further 
        np.savez(self.get_output('gaucov_sims')[:-4],cl_sims=cells_sims)
        
        #Compute covariance
        covar=np.cov(cells_sims.T)
        return covar

    def get_covar_analytic(self,lth,clth,bpws,wsp) :
        #Create a dummy file for the covariance MCM
        f=open(self.get_output('gaucov_sims'),"w")
        f.close()

        covar=np.zeros([self.ncross*self.nell,self.ncross*self.nell])
        cwsp=self.get_covar_mcm(wsp)

        ix_1=0
        for i1 in range(self.nbins) :
            for j1 in range(i1,self.nbins) :
                ix_2=0
                for i2 in range(self.nbins) :
                    for j2 in range(i2,self.nbins) :
                        ca1b1=clth[self.ordering[i1,i2]]
                        ca1b2=clth[self.ordering[i1,j2]]
                        ca2b1=clth[self.ordering[j1,i2]]
                        ca2b2=clth[self.ordering[j1,j2]]
                        cov_here=nmt.gaussian_covariance_flat(cwsp,lth,ca1b1,ca1b2,ca2b1,ca2b2)
                        covar[ix_1*self.nell:(ix_1+1)*self.nell,:][:,ix_2*self.nell:(ix_2+1)*self.nell]=cov_here
                        ix_2+=1
                ix_1+=1

        return covar

    def get_masks(self) :
        #Depth-based mask
        self.fsk,mp_depth=read_flat_map(self.get_input("depth_map"),i_map=0)
        mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
        msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=self.config['depth_cut']]=1

        fskb,mskfrac=read_flat_map(self.get_input("masked_fraction"),i_map=0)
        compare_infos(self.fsk,fskb)
        
        #Create binary mask (fraction>threshold and depth req.)
        msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>self.config['mask_thr']]=1
        msk_bi=msk_bo*msk_depth

        if self.config['mask_systematics'] :
            '''
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
            '''
            raise NotImplementedError("TODO: implement systematics masking")
        return msk_bi,mskfrac,mp_depth

    def get_tracers(self,temps) :
        hdul=fits.open(self.get_input('ngal_maps'))
        if len(hdul)%2!=0 :
            raise ValueError("Input file should have two HDUs per map")
        nbins=len(hdul)//2
        tracers_nocont=[Tracer(hdul,i,self.fsk,self.msk_bi,self.mskfrac,contaminants=None)
                        for i in range(nbins)]
        tracers_wcont=[Tracer(hdul,i,self.fsk,self.msk_bi,self.mskfrac,contaminants=temps)
                       for i in range(nbins)]
        hdul.close()
        return tracers_nocont,tracers_wcont

    def get_contaminants(self) :
        temps=[]
        #Depth
        temps.append(self.mp_depth)
        #Dust
        for t in self.read_map_bands(self.get_input('dust_map'),False,self.config['band']) :
            temps.append(t)
        #Stars
        fskb,t=read_flat_map(self.get_input('star_map'),i_map=0)
        compare_infos(self.fsk,fskb)
        temps.append(t)
        #Observing conditions
        #TODO: we were not marginalizing over all of these
        #for oc in ['airmass','ccdtemp','ellipt','exptime','nvisit',
        #           'seeing','sigma_sky','skylevel'] :
        for oc in ['airmass','seeing','sigma_sky'] :
            for t in self.read_map_bands(self.get_input(oc+'_maps'),
                                         self.config['oc_all_bands'],self.config['band']) :
                temps.append(t)
        temps=np.array(temps)
        #Remove mean
        for i_t,t in enumerate(temps) :
            temps[i_t]-=np.sum(self.msk_bi*self.mskfrac*t)/np.sum(self.msk_bi*self.mskfrac)

        return temps

    def parse_input(self) :
        #TODO: Parse input params

        return

    def get_sacc_tracers(self,tracers) :
        sacc_tracers=[]
        for i_t,t in enumerate(tracers) :
            z=(t.nz_data['z_i']+t.nz_data['z_f'])*0.5
            nz=t.nz_data['nz_cosmos']
            T=sacc.Tracer('bin_%d'%i_t,'point',z,nz,exp_sample="HSC_DESC")
            T.addColumns({'nz_'+c:t.nz_data['nz_'+c] 
                          for c in ['demp','ephor','ephor_ab','frankenz','nnpz']})
            sacc_tracers.append(T)

        return sacc_tracers

    def write_vector_to_sacc(self,fname_out,sacc_t,sacc_b,cls,covar=None,verbose=False) :
        sacc_mean=sacc.MeanVec(cls.flatten())
        if covar is None :
            sacc_precision=None
        else :
            sacc_precision=sacc.Precision(covar,"dense",is_covariance=True, binning=sacc_b)

        sacc_meta={'Area_rad':self.area_patch}
        s=sacc.SACC(sacc_t,sacc_b,sacc_mean,precision=sacc_precision,meta=sacc_meta)
        if verbose :
            s.printInfo()
        s.saveToHDF(fname_out)

    def get_sacc_binning(self,ell_eff,lini,lend,windows=None) :
        typ,ell,dell,t1,q1,t2,q2=[],[],[],[],[],[],[]
        for t1i in range(self.nbins) :
            for t2i in range(t1i,self.nbins) :
                for i_l,l in enumerate(ell_eff) :
                    typ.append('F')
                
                    ell.append(l)
                    dell.append(lend[i_l]-lini[i_l])
                    t1.append(t1i)
                    q1.append('P')
                    t2.append(t2i)
                    q2.append('P')

        if windows is None :
            return sacc.Binning(typ,ell,t1,q1,t2,q2,deltaLS=dell)
        else :
            return sacc.Binning(typ,ell,t1,q1,t2,q2,deltaLS=dell,windows=windows)
                
    def run(self) :
        self.parse_input()

        print("Reading mask")
        self.msk_bi,self.mskfrac,self.mp_depth=self.get_masks()

        print("Computing area")
        self.area_pix=np.radians(self.fsk.dx)*np.radians(self.fsk.dy)
        self.area_patch=np.sum(self.msk_bi*self.mskfrac)*self.area_pix
        self.lmax=int(180.*np.sqrt(1./self.fsk.dx**2+1./self.fsk.dy**2))

        print("Reading contaminants")
        temps=self.get_contaminants()

        print("Setting bandpowers")
        lini=np.array(self.config['ell_bpws'])[:-1]
        lend=np.array(self.config['ell_bpws'])[ 1:]
        bpws=nmt.NmtBinFlat(lini,lend)
        ell_eff=bpws.get_effective_ells()

        print("Generating tracers")
        tracers_nc,tracers_wc=self.get_tracers(temps)
        self.nbins=len(tracers_nc)

        print("Translating into SACC tracers")
        tracers_sacc=self.get_sacc_tracers(tracers_nc)

        self.ordering=np.zeros([self.nbins,self.nbins],dtype=int)
        ix=0
        for i in range(self.nbins) :
            for j in range(i,self.nbins) :
                self.ordering[i,j]=ix
                if j!=i :
                    self.ordering[j,i]=ix

        print("Getting MCM")
        wsp=self.get_mcm(tracers_nc,bpws)

        print("Computing window function")
        windows=self.get_sacc_windows(wsp)

        print("Computing SACC binning")
        #No windows
        binning_nw=self.get_sacc_binning(ell_eff,lini,lend,windows=None)
        #With windows
        binning_ww=self.get_sacc_binning(ell_eff,lini,lend,windows=windows)

        print("Computing power spectra")
        cls_wodpj,_=self.get_power_spectra(tracers_nc,wsp,bpws)
        cls_wdpj,cls_wdpj_coupled=self.get_power_spectra(tracers_wc,wsp,bpws)
        self.ncross,self.nell=cls_wodpj.shape

        print("Getting guess power spectra")
        lth,clth=self.get_cl_guess(ell_eff,cls_wdpj)

        print("Computing deprojection bias")
        cls_wdpj,cls_deproj=self.get_dpj_bias(tracers_wc,lth,clth,cls_wdpj_coupled,wsp,bpws)

        print("Computing covariance")
        cov_wodpj=self.get_covar(lth,clth,bpws,wsp,None,None)
        if self.config['gaus_covar_type']=='analytic' :
            cov_wdpj=cov_wodpj
        else :
            cov_wdpj=self.get_covar(lth,clth,bpws,wsp,temps,cls_deproj)

        if self.config['add_ssc'] :
            raise NotImplementedError("SSC not working yet")
            print("Computing SSC")
            ng_data,cov_ssc=self.get_covar_ssc(tracers_nc,ell_eff)
            cov_wodpj+=cov_ssc
            cov_wdpj+=cov_ssc
            print(ng_data)
            #np.savetxt(o.prefix_out+".ngals",ng_data,
            #           header='[1]-Ngals [2]-Sigma_poisson [3]-Sigma_CV [4]-Sigma_T')

        print("Computing noise bias")
        nls=self.get_noise(tracers_nc,wsp,bpws)

        print("Writing output")
        self.write_vector_to_sacc(self.get_output('noi_bias'),tracers_sacc,
                                  binning_nw,nls,verbose=False)
        self.write_vector_to_sacc(self.get_output('dpj_bias'),tracers_sacc,
                                  binning_nw,cls_deproj,verbose=False)
        self.write_vector_to_sacc(self.get_output('power_spectra_wodpj'),tracers_sacc,
                                  binning_ww,cls_wodpj,covar=cov_wodpj,verbose=False)
        self.write_vector_to_sacc(self.get_output('power_spectra_wdpj'),tracers_sacc,
                                  binning_ww,cls_wdpj,covar=cov_wdpj,verbose=True)
        
        exit(1)

if __name__ == '__main__':
    cls = PipelineStage.main()
