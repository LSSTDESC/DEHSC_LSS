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

prefix_data='/global/cscratch1/sd/damonge/HSC/'
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
parser.add_option('--output-file', dest='fname_out',default=None,type=str,
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

####
# Read options
(o, args) = parser.parse_args()

print("Reading mask")
#Create depth-based mask
fsk,mp_depth=fm.read_flat_map(o.prefix_in+"_10s_depth_mean_fluxerr.fits",2)
mp_depth[np.isnan(mp_depth)]=0; mp_depth[mp_depth>40]=0
msk_depth=np.zeros_like(mp_depth); msk_depth[mp_depth>=o.depth_cut]=1

#Read masked fraction
fskb,mskfrac=fm.read_flat_map(o.prefix_in+'_MaskedFraction.fits',i_map=0)
fm.compare_infos(fsk,fskb)

#Create BO-based mask
msk_bo=np.zeros_like(mskfrac); msk_bo[mskfrac>o.mask_thr]=1
msk_t=msk_bo*msk_depth

#Read contaminants
print("Reading contaminant templates")
temps=[]
def read_map_bands(fname,read_bands) :
  if read_bands :
    i_map=-1
  else :
    i_map=['g','r','i','z','y'].index(o.band)
  fskb,temp=fm.read_flat_map(fname,i_map=i_map)
  fm.compare_infos(fsk,fskb)
  if i_map!=-1 :
    temp=[temp]

  return temp
# 1- Depth
if o.cont_depth :
  temps.append(mp_depth)
# 2- Dust
if o.cont_dust :
  temp=read_map_bands(o.prefix_in+'_syst_dust.fits',o.cont_dust_bands)
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
    temp=read_map_bands(o.prefix_in+'_oc_'+c+'.fits',o.cont_oc_bands)
    for t in temp :
      temps.append(t)
if temps==[] :
  temps=None
else :
  print(" - Will marginalize over a total of %d contaminant templates"%(len(temps)))

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
    self.weight=masked_fraction*mask_binary
    goodpix=np.where(mask_binary>0.1)[0]
    ndens=np.sum(nmap*mask_binary)/np.sum(self.weight)
    self.delta=np.zeros_like(self.weight)
    self.delta[goodpix]=nmap[goodpix]/(ndens*masked_fraction[goodpix])-1

    #Reshape contaminants
    conts=None
    if contaminants is not None :
      conts=[[c.reshape([self.fsk.ny,self.fsk.nx]) for c in contaminants]]

    #Form NaMaster field
    self.field=nmt.NmtFieldFlat(np.radians(self.fsk.lx),np.radians(self.fsk.ly),
                                self.weight.reshape([self.fsk.ny,self.fsk.nx]),
                                [self.delta.reshape([self.fsk.ny,self.fsk.nx])],
                                templates=conts)

hdul=fits.open(o.fname_maps)
if len(hdul)%2!=0 :
  raise ValueError("Input file should have two HDUs per map")
nbins=len(hdul)/2
tracers=[Tracer(hdul,i,fsk,msk_t,mskfrac,contaminants=temps) for i in np.arange(nbins)]
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

#Compute all cross-correlations
print("Computing all cross-power specra")
print("Warning : deprojection bias still missing")
cls_all=[]
for i in range(nbins) :
  for j in range(i,nbins) :
    cl_coupled=nmt.compute_coupled_cell_flat(tracers[i].field,tracers[j].field,bpws)
    cls_all.append(wsp.decouple_cell(cl_coupled)[0])
cls_all=np.array(cls_all)

#Save to SACC format
print("Saving to SACC")
#Tracers
sacc_tracers=[]
for i_t,t in enumerate(tracers) :
  z=(t.nz_data['z_i']+t.nz_data['z_f'])*0.5
  nz=t.nz_data['n_z']
  T=sacc.Tracer('bin_%d'%i_t,'point',z,nz,exp_sample=o.hsc_field)
  sacc_tracers.append(T)
#Binning and mean
type,ell,dell,t1,q1,t2,q2=[],[],[],[],[],[],[]
i_cross=0
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
sacc_binning=sacc.Binning(type,ell,t1,q1,t2,q2,deltaLS=dell)
sacc_mean=sacc.MeanVec(cls_all.flatten())
s=sacc.SACC(sacc_tracers,sacc_binning,sacc_mean)
s.printInfo()
s.saveToHDF(o.fname_out)
