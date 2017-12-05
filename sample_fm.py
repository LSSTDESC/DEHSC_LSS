import numpy as np
import matplotlib.pyplot as plt
import flatmaps as fm
import pandas as pd
from astropy.io import fits
from astropy.table import Table

fname="/global/cscratch1/sd/damonge/HSC/HSC_WIDE_WIDE_GAMA15H_forced.fits"

#Create information for a flat-sky map covering the GAMA15H field with resolution ~3.4 arcmin
mi=fm.FlatMapInfo([212.5,222.],[-2.,2.],dx=0.057,dy=0.057)

print "Reading"
data=(Table.read(fname,format='fits')).to_pandas()

print "Cleaning" #Basically a copy from Humna's cleaning steps
isNulls= [k for k in data.keys() if k.__contains__('isnull')] # all the columns that signify nulls
indToDrop= []
for col in isNulls:
    nullInd= np.where(data[col]==True)[0]
    indToDrop+=nullInd
# find unique indexes
indToDrop= list(set(indToDrop))
if len(indToDrop)>0: # remove the objs. with any null entries
    print 'Dropped %s entries.'%len(indToDrop)
    data= data.drop(indToDrop)
# drop the columns now. irrelevant.
data= data.drop(isNulls, axis= 1)
print 'Dropped %s columns.'%len(isNulls)
noClassInd= np.where(np.isnan(data['iclassification_extendedness']))[0]
print 'Dropped %s entries.'%len(noClassInd)
data= data.drop(noClassInd, axis= 0)
print 'Final size: ', np.shape(data)
before= len(data)
data= data.dropna(axis= 0)
print 'Dropped %s rows'%(before-len(data.keys()))
print 'Final size: ', np.shape(data)

#Make a map from galaxy positions
print "Pixel indices"
ipix=mi.pos2pix(data['ra'],data['dec'])
print "Binning"
#Bin positions into pixels
mp=np.bincount(ipix,minlength=mi.get_size())

print "Writing map"
mi.write_flat_map("map_test",mp)

print "Reading map"
mi2,mp2=fm.read_flat_map("map_test.npz")

print "Plotting"
#Plot resulting map
mi.view_map(mp)
mi2.view_map(mp2)
plt.show()
